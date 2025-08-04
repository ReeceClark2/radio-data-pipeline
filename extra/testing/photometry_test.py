from pathlib import Path
from rc_test import run_cartographer


def extract_corrected_photometry(output_txt_path: Path) -> float:
    with open(output_txt_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    data_lines = [line for line in lines if line.startswith("RAPIX") or line[0].isdigit()]
    if len(data_lines) < 2:
        raise ValueError(f"{output_txt_path} does not contain header and data rows.")

    header = data_lines[0].split()
    data = data_lines[1].split()

    try:
        idx = header.index("Correctd_Phtmtry")
        return float(data[idx])
    except (ValueError, IndexError) as e:
        raise ValueError(f"Could not extract Correctd_Phtmtry from {output_txt_path}") from e


def scan_photometry_results(root_dir: Path) -> dict:
    results = {}
    for txt_path in root_dir.rglob("Output.txt"):
        try:
            folder_name = txt_path.parent.name  # e.g. "0137060_phot"
            base_name = folder_name.replace("_phot", "")
            value = extract_corrected_photometry(txt_path)
            results[base_name] = value
        except Exception as e:
            print(f"Warning: {e}")
    return results



# Example usage
if __name__ == "__main__":
    results_dir = Path("results")
    values = scan_photometry_results(results_dir)
    for name, value in values.items():
        print(f"{name}: {value}")
