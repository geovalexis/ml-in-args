import pandas as pd

if __name__ == "__main__":
    orig_batch = pd.read_csv("results/args_table_batch2.csv", header=0)
    dest_batch = pd.read_csv("results/args_table_only_batch3.csv", header=0)
    all_df = pd.concat([orig_batch, dest_batch], ignore_index=True)
    all_df.fillna(0, inplace=True)
    all_df.to_csv("results/args_table_batch3 .csv", index=False)