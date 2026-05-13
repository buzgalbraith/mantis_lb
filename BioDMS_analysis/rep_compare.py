import numpy as np

data = {
        "Files_RR":   [150, 150, 150, 150, 150],
            "Files_Sim":  [178, 149, 151, 132, 140],
            "Files_Asym": [180, 142, 152, 133, 144],
            "Kmers_RR":   [10.20, 7.46, 7.19, 6.96, 6.96],  # x10^7
            "Kmers_Sim":  [9.11,  6.61, 6.72, 6.89, 6.76],
            "Kmers_Asym": [10.80, 7.32, 7.56, 6.30, 6.52],
            "EC_RR":   [6.46, 6.39, 6.50, 6.34, 6.05],       # x10^6
            "EC_Sim":  [5.37, 6.34, 6.57, 5.54, 5.86],
            "EC_Asym": [8.32, 5.90, 6.24, 5.43, 5.51],
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
        "Files":          ("Files_RR", "Files_Sim", "Files_Asym"),
            "Kmers (x10^7)":  ("Kmers_RR", "Kmers_Sim", "Kmers_Asym"),
            "Equiv. Classes (x10^6)": ("EC_RR", "EC_Sim", "EC_Asym"),
}

for group_name, (rr_key, sim_key, asym_key) in groups.items():
        ## total ## 
        rr_total   = sum(data[rr_key])
        sim_total  = sum(data[sim_key])
        asym_total = sum(data[asym_key])

        sim_pct  = (sim_total  - rr_total) / rr_total * 100
        asym_pct = (asym_total - rr_total) / rr_total * 100
        ## mean ## 
        rr_mean = np.mean(data[rr_key])
        sim_mean  = np.mean(data[sim_key])
        asym_mean = np.mean(data[asym_key])

        sim_pct_mean  = (sim_mean  - rr_mean) / rr_mean * 100
        asym_pct_mean = (asym_mean - rr_mean) / rr_mean * 100
        ## variance ## 
        rr_variance = np.var(data[rr_key])
        sim_variance  = np.var(data[sim_key])
        asym_variance = np.var(data[asym_key])

        sim_pct_variance  = (sim_variance  - rr_variance) / rr_variance * 100
        asym_pct_variance = (asym_variance - rr_variance) / rr_variance * 100

        print(f"\n{group_name}")
        print(f"  RR total:   {rr_total:.2f}")
        print(f"  Sim total:  {sim_total:.2f}  ({sim_pct:+.2f}% vs RR)")
        print(f"  Asym total: {asym_total:.2f}  ({asym_pct:+.2f}% vs RR)")
        print(f"  RR mean:   {rr_mean:.2f}")
        print(f"  Sim mean:  {sim_mean:.2f}  ({sim_pct_mean:+.2f}% vs RR)")
        print(f"  Asym mean: {asym_mean:.2f}  ({asym_pct_mean:+.2f}% vs RR)")
        print(f"  RR variance:   {rr_variance:.2f}")
        print(f"  Sim variance:  {sim_variance:.2f}  ({sim_pct_variance:+.2f}% vs RR)")
        print(f"  Asym variance: {asym_variance:.2f}  ({asym_pct_variance:+.2f}% vs RR)")
