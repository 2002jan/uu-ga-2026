import os
import csv
from pathlib import Path

results_path = Path("./results")
data_root = results_path / "htseq_counts"
output_root = results_path / "htseq_combined"

output_root.mkdir(exist_ok=True)

conditions = ["BH", "Serum"]

for c in conditions:
    c_path = data_root / c

    c_samples = set()

    for f in c_path.iterdir():
        c_samples.add(f.name.split(".")[0])  

    for s in c_samples:
        p_path = c_path / f"{s}.counts_paired.tsv"
        s_path = c_path / f"{s}.counts_single.tsv"

        if not s_path.is_file() or not p_path.is_file():
            continue

        combined = []

        with open(p_path, 'r') as p_file:
            with open(s_path, 'r') as s_file:
                p_tsv = csv.reader(p_file, delimiter="\t")
                s_tsv = csv.reader(s_file, delimiter="\t")

                for p_r, s_r in zip(p_tsv, s_tsv):
                    if(p_r[0] != s_r[0]):
                        print("No match")
                        exit(1)

                    combined.append([p_r[0], p_r[1], int(p_r[2]) + int(s_r[2])])

        combined_path = output_root / f"{s}_{c}.counts.tsv"

        if combined_path.exists():
            print(f"{combined_path} already exists... deleting...")
            combined_path.unlink()

        with open(combined_path, 'w') as c_file:
            writer = csv.writer(c_file, delimiter="\t")
            writer.writerows(combined)

combined_path = output_root / f"combined.tsv"

if combined_path.exists():
    print("Output file exists! Deleting")
    combined_path.unlink()

combined_rows = []
header = ["Gene Id", "Name"]

for f in output_root.iterdir():
    if not f.name.endswith(".counts.tsv"):
        continue
    
    name = f.name.split(".")[0]
    header.append(name)

    with open(f, 'r') as tsv_f:
        reader = csv.reader(tsv_f, delimiter="\t")
        i = 0
        firstRun = len(combined_rows) < 1

        for row in reader:
            if firstRun:
                combined_rows.append(row)
            else:
                current_row = combined_rows[i]

                if current_row[0] != row[0]:
                    print("Id mismatch")
                    print(current_row[0], row[0])
                    exit(1)

                combined_rows[i].append(row[2])
            i += 1

with open(combined_path, 'w') as output_file:
    writer = csv.writer(output_file, delimiter="\t")

    writer.writerow(header)
    writer.writerows(combined_rows)

print(f"Done... results saved to: {combined_path}")
