import numpy as np

data = {
        "CQF_RR":   [2.85, 2.85, 2.85, 2.85, 2.85],
            "CQF_Sim":  [2.85, 2.85, 2.85, 2.85, 2.85],
            "CQF_Asym": [2.85, 2.85, 2.85, 2.85, 2.85],
            "CDBG_RR":   [79.5, 79.5, 83.0, 80.1, 75.7],
            "CDBG_Sim":  [71.2, 81.0, 87.2, 65.3, 69.3],
            "CDBG_Asym": [119.2, 70.3, 79.3, 62.1, 66.6],
            "MST_RR":   [29.3, 31.0, 31.9, 30.8, 29.2],
            "MST_Sim":  [26.4, 31.6, 32.9, 26.8, 28.9],
            "MST_Asym": [41.3, 28.6, 30.9, 26.1, 26.5],
}

print(f"{'Column':<12} {'Mean':>10} {'Variance':>12} {'Std Dev':>10}")
print("-" * 46)
for col, vals in data.items():
        arr = np.array(vals)
        print(f"{col:<12} {arr.mean():>10.4f} {arr.var():>12.4f} {arr.std():>10.4f}")

    # % change from RR for each metric group
print("\n% Change from Round Robin (by total)")
print("-" * 40)
groups = {
        "CQF (GB)":  ("CQF_RR",  "CQF_Sim",  "CQF_Asym"),
            "CDBG (MB)": ("CDBG_RR", "CDBG_Sim", "CDBG_Asym"),
            "MST (MB)":  ("MST_RR",  "MST_Sim",  "MST_Asym"),
}

for group_name, (rr_key, sim_key, asym_key) in groups.items():
        rr_total   = sum(data[rr_key])
        sim_total  = sum(data[sim_key])
        asym_total = sum(data[asym_key])

        sim_pct  = (sim_total  - rr_total) / rr_total * 100
        asym_pct = (asym_total - rr_total) / rr_total * 100

        print(f"\n{group_name}")
        print(f"  RR total:   {rr_total:.2f}")
        print(f"  Sim total:  {sim_total:.2f}  ({sim_pct:+.2f}% vs RR)")
        print(f"  Asym total: {asym_total:.2f}  ({asym_pct:+.2f}% vs RR)")
