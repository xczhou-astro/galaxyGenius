import json
import os
import sys
import time

import numpy as np
import requests
from astropy.cosmology import Planck15


TNG_API_BASE = "https://www.tng-project.org/api/"


def make_request(url, headers, max_retries=5):
    start_time = time.time()
    content = None

    for attempt in range(max_retries):
        try:
            response = requests.get(url, headers=headers, timeout=120)
            response.raise_for_status()
            content = response
            break
        except requests.exceptions.RequestException as e:
            if attempt == max_retries - 1:
                print(f"Error: Failed to make request after {max_retries} attempts.")
                print(f"URL: {url}")
                print(f"Error message: {str(e)}")
            else:
                print(
                    f"Request failed (attempt {attempt + 1}/{max_retries}). Retrying..."
                )
                time.sleep(2**attempt)

    print(f"Requesting time taken: {time.time() - start_time:.2f} seconds")

    if content is None:
        sys.exit(1)

    return content


def get_snap_info(snapshot_id, stellarMass_ranges, headers, cosmo, tng_version, save_path):
    h = cosmo.h
    min_stellar_mass = np.float32(stellarMass_ranges[0]) / 10**10 * h

    url = (
        f"{TNG_API_BASE}{tng_version}/snapshots/{snapshot_id}"
        f"/subhalos/?mass_stars__gt={min_stellar_mass}&subhaloflag=1&limit=1000000"
    )
    if np.isfinite(stellarMass_ranges[1]):
        max_stellar_mass = np.float32(stellarMass_ranges[1]) / 10**10 * h
        url += f"&mass_stars__lt={max_stellar_mass}"

    print("snapshot url: ", url)
    response = make_request(url, headers=headers)
    data = response.json()
    results = data.get("results", [])

    print("Number of subhalos: ", len(results))

    with open(os.path.join(save_path, f"snapshot-{snapshot_id}.json"), "w") as f:
        json.dump(results, f, indent=4)

    return results


def _subhalo_detail_url(subhalo_id, snapshot_id, tng_version):
    return (
        f"{TNG_API_BASE}{tng_version}/snapshots/{snapshot_id}/subhalos/{subhalo_id}/"
    )


def get_subhalo_info(subhalo_id, headers, snapshot_id, tng_version, results=None):
    """
    Fetch subhalo JSON. If ``results`` is given and contains ``subhalo_id``,
    reuse the catalog URL; otherwise request the canonical detail URL (works
    for arbitrary IDs not present in a mass-filtered listing).
    """
    subhalo_url = None
    if results:
        id_to_result = {r["id"]: r for r in results}
        if subhalo_id in id_to_result:
            subhalo_url = id_to_result[subhalo_id]["url"]

    if subhalo_url is None:
        subhalo_url = _subhalo_detail_url(subhalo_id, snapshot_id, tng_version)

    print("subhalo url: ", subhalo_url)
    response = make_request(subhalo_url, headers=headers)
    return response.json()


def request_subhalo_particles(cutout_url, headers, savefilename, max_retries=5):
    for attempt in range(max_retries):
        try:
            with requests.get(
                cutout_url, headers=headers, stream=True, timeout=300
            ) as response:
                response.raise_for_status()
                with open(savefilename, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            return
        except requests.exceptions.RequestException:
            print(f"Attempt {attempt + 1}/{max_retries} failed, retrying...")
            if os.path.isfile(savefilename):
                try:
                    os.remove(savefilename)
                except OSError:
                    pass
            if attempt < max_retries - 1:
                time.sleep(2**attempt)
            else:
                raise


def get_particles(subhalo_id, snapshot_id, headers, save_path, tng_version):
    cutout_url = (
        f"{TNG_API_BASE}{tng_version}/snapshots/{snapshot_id}/subhalos/"
        f"{subhalo_id}/cutout.hdf5"
    )
    print("cutout url: ", cutout_url)
    out_path = os.path.join(save_path, f"subhalo_{subhalo_id}_particles.h5")
    request_subhalo_particles(cutout_url, headers, out_path)


if __name__ == "__main__":
    snapshot_ids = [94]
    snapshot_zs = [0.06]
    subhalo_ids = [0, 31, 253881]
    # subhalo_ids = [31]
    min_stellar_mass = 10**10
    max_stellar_mass = np.inf
    tng_version = "TNG100-1"
    save_path_base = 'TNGData'

    stellarMass_ranges = [min_stellar_mass, max_stellar_mass]

    cosmo = Planck15
    api_key = os.environ.get("TNG_API_KEY")
    if not api_key:
        print("Set TNG_API_KEY in the environment (IllustrisTNG API token).", file=sys.stderr)
        sys.exit(1)
    headers = {"api-key": api_key}

    with open("failed_subhalos.txt", "w") as f:
        f.write("snapshot-ID,subhaloID\n")

    for snap_id, snap_z in zip(snapshot_ids, snapshot_zs):
        print("Processing snapshot", snap_id, "at redshift", snap_z)

        save_path = f"{save_path_base}/{tng_version}/snapshots-{snap_id}"
        os.makedirs(save_path, exist_ok=True)

        results = get_snap_info(
            snap_id, stellarMass_ranges, headers, cosmo, tng_version, save_path
        )

        ids_to_process = (
            subhalo_ids if subhalo_ids is not None else [r["id"] for r in results]
        )

        for subhalo_id in ids_to_process:
            sub_dir = f"{save_path_base}/{tng_version}/snapshots-{snap_id}/subhalo-{subhalo_id}"
            os.makedirs(sub_dir, exist_ok=True)

            particle_path = os.path.join(
                sub_dir, f"subhalo_{subhalo_id}_particles.h5"
            )
            if os.path.exists(particle_path):
                print(f"Subhalo {subhalo_id} already downloaded")
                continue

            data = get_subhalo_info(
                subhalo_id, headers, snap_id, tng_version, results=results
            )

            with open(os.path.join(sub_dir, f"subhalo_{subhalo_id}.json"), "w") as f:
                json.dump(data, f, indent=4)

            try:
                get_particles(subhalo_id, snap_id, headers, sub_dir, tng_version)
            except Exception as e:
                print("Error getting particles: ", e)
                with open("failed_subhalos.txt", "a") as f:
                    f.write(f"{snap_id},{subhalo_id}\n")
