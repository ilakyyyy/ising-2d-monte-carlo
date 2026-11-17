import glob
import re

output_file = "results_L_16_128/combined_simulations_L16to128.txt"
files = glob.glob("results_L_16_128/simulations_L*.txt")

with open(output_file, 'w') as outfile:
    outfile.write(f"{'L':<8} {'beta':<10} {'energy':<12} {'spec-heat':<12} {'abs-magnet':<12} {'suscept'}\n")

    for filename in files:
        match = re.search(r'simulations_L([\d.]+)\.txt', filename)
        if not match:
            continue

        l_value = match.group(1)
        with open(filename, 'r') as infile:
            for line in infile.readlines()[1:]:
                line = line.strip()
                if line:
                    outfile.write(f"{l_value:<8} {line}\n")

print(f"Combined {len(files)} files into {output_file}")
