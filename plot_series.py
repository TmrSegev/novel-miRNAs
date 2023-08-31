import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_series(series, ticks):
    series.plot.hist()
    plt.title(series.name)
    print("name:", series.name, "min:", series.min(), "max:", series.max())
    plt.xticks(np.arange(series.min(), series.max() + ticks, ticks), rotation=45)
    mean = series.mean()
    std = series.std()
    plt.axvline(mean, color="red")
    plt.axvline(mean + std, color="yellow")
    plt.axvline(mean - std, color="yellow")
    plt.savefig("./figures/mirgenedb_{}.png".format(series.name), dpi=300)
    plt.clf()

if __name__ == '__main__':
    mirgenedb = pd.read_csv("temp.csv", sep='\t')
    plot_series(mirgenedb['Hairpin_seq_trimmed_length'], 5.0)
    plot_series(mirgenedb['Mature_connections'], 1.0)
    plot_series(mirgenedb['Mature_BP_ratio'].astype('float'), 0.05)
    plot_series(mirgenedb['Mature_max_bulge'].astype('float'), 1.0)
    plot_series(mirgenedb['Loop_length'], 2.0)
    plot_series(mirgenedb['Mature_Length'], 1.0)
    plot_series(mirgenedb['Star_length'], 1.0)
    plot_series(mirgenedb['Star_connections'], 1.0)
    plot_series(mirgenedb['Star_BP_ratio'].astype('float'), 0.05)
    plot_series(mirgenedb['Star_max_bulge'].astype('float'), 1.0)
    plot_series(mirgenedb['Max_bulge_symmetry'], 1.0)
    plot_series(mirgenedb['min_one_mer_hairpin'], 0.05)
    plot_series(mirgenedb['max_one_mer_hairpin'], 0.05)

