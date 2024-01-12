import pickle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

sdss_gor = pickle.load(open('sdss_gor.p', 'rb'))

sdss_sample_medians = []

for sample in sdss_gor:
    this_sample = sdss_gor[sample]

    all_gor = []
    for r_vals, this_gor in this_sample:
        all_gor.append(this_gor)
        plt.plot(r_vals, this_gor, alpha=0.2, color='red')

    all_gor = np.array(all_gor)

    gor_med = np.median(all_gor, axis=0)

    sdss_sample_medians.append(gor_med)

    plt.plot(r_vals, gor_med, color='red', lw=2, label='median')
    plt.legend()

    plt.xlabel(r'$r$ ($h^{-1}$ Mpc)', fontsize=16)
    plt.ylabel(r'$g(r)$', fontsize=16)

    plt.savefig(f'figures/SDSS_{sample}_gor.jpg', dpi=600)
    plt.cla()
    plt.clf()


# Plot medians for ONLY SDSS sample1
plt.plot(r_vals, sdss_sample_medians[1], color='red', zorder=1000, lw=3)

mocks_gor = pickle.load(open('mocks_gor.p', 'rb'))

mocks_medians = []
gor_med = []

for mock in mocks_gor:
    this_mock = mocks_gor[mock]

    all_gor = []
    for r_vals, this_gor in this_mock:
        all_gor.append(this_gor)

    all_gor = np.array(all_gor)

    gor_med = np.median(all_gor, axis=0)

    plt.plot(r_vals, gor_med, color='blue', alpha=0.2)


legend_elements = [Line2D([0], [0], color='r', lw=3, label='SDSS median')\
                , Line2D([0], [0], color='b', lw=1, label='mocks medians')]

plt.legend(handles=legend_elements)

plt.xlabel(r'$r$ ($h^{-1}$ Mpc)', fontsize=16)
plt.ylabel(r'$g(r)$', fontsize=16)

plt.tight_layout()

plt.savefig('figures/SDSS_mocks_gor_medians.jpg', dpi=600)


