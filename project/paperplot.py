from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import matplotlib.ticker as mticker


def read_depth_file(addpath):
    folder_path = os.getcwd() + addpath
    xvalue = []
    yvalue = []
    for filename in os.listdir(folder_path):
        # print(filename)
        if filename.startswith("W=15_D="):
            file_path = os.path.join(folder_path, filename)
            d_value = float(filename.split('_')[1][2:])
            with open(file_path, 'r') as file:
                file_content = file.read().strip()
            xvalue.append(d_value)
            yvalue.append(file_content)

    zipped_pairs = zip(xvalue, yvalue)
    sorted_pairs = sorted(zipped_pairs)

    xvalue, yvalue = zip(*sorted_pairs)
    xvalue = np.array(xvalue).astype(float)
    yvalue = np.array(yvalue).astype(float)
    xvalue *= np.log(xvalue) * np.log(xvalue)
    return xvalue, yvalue


def read_width_file(addpath, number):
    folder_path = os.getcwd() + addpath
    xvalue = []
    yvalue = []
    if number == 20:
        pattern = r"W=.*_D=20_.*"
    elif number == 30:
        pattern = r"W=.*_D=30_.*"
    for filename in os.listdir(folder_path):
        if re.match(pattern, filename):
            print(filename)
            file_path = os.path.join(folder_path, filename)
            w_value = float(filename.split('_')[0][2:])
            with open(file_path, 'r') as file:
                file_content = file.read().strip()
            xvalue.append(w_value)
            yvalue.append(file_content)

    zipped_pairs = zip(xvalue, yvalue)
    sorted_pairs = sorted(zipped_pairs)

    xvalue, yvalue = zip(*sorted_pairs)
    xvalue = np.array(xvalue).astype(float)
    yvalue = np.array(yvalue).astype(float)
    # xvalue *= np.log(xvalue) * np.log(xvalue)
    return xvalue, yvalue


def draw1():

    xvalue_gradient, yvalue_gradient = read_depth_file('/data_toymodel/')
    xvalue_pulse, yvalue_pulse = read_depth_file('/data_toymodel_pulse/')

    coefficients = np.polyfit(
        xvalue_gradient, yvalue_gradient, 1)  # 1 表示一次多项式（线性）
    slope1, intercept1 = coefficients
    y_fit_gradient = slope1 * xvalue_pulse + intercept1

    coefficients = np.polyfit(xvalue_pulse, yvalue_pulse, 1)  # 1 表示一次多项式（线性）
    slope2, intercept2 = coefficients
    y_fit_pulse = slope2 * xvalue_pulse + intercept2

    # fig = plt.figure(dpi=200, figsize=(20, 5.5))
    fig = plt.figure(dpi=200)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'Times New Roman'

    ax1 = plt.subplot(1, 1, 1)
    ax1.scatter(xvalue_gradient, yvalue_gradient, marker='o', edgecolors='#9a031e',
                facecolors='none', linewidth=3, s=60, label='Energy Gradient', zorder=5)
    ax1.plot(xvalue_pulse, y_fit_gradient, '--', c='#03045e',
             linewidth=2.5, zorder=5, label='Fit of Energy Gradient')
    ax1.scatter(xvalue_pulse, yvalue_pulse, marker='o', edgecolors='#fcbf49',
                facecolors='none', linewidth=3, s=60, label='Wave Packet')
    ax1.plot(xvalue_pulse, y_fit_pulse, '-', c='#0096c7',
             linewidth=2.5, zorder=3, label='Fit of Wave Packet')

    ax = plt.gca()
    mylinewidth = 1.5
    ax.spines['bottom'].set_linewidth(mylinewidth)
    ax.spines['left'].set_linewidth(mylinewidth)
    ax.spines['top'].set_linewidth(mylinewidth)
    ax.spines['right'].set_linewidth(mylinewidth)

    plt.xlim(left=0, right=550)
    plt.ylim(bottom=0, top=65000)

    yticks_name = ['0'] + [f'{y}'+r'$\times 10^4$' for y in [
        1, 2, 3, 4, 5, 6]]
    yticks = [y*1e4 for y in [0, 1, 2, 3, 4, 5, 6]]
    # print(yticks_name)
    plt.xticks(fontsize=12)
    plt.yticks(yticks, yticks_name, fontsize=12, fontname='Times New Roman')

    plt.xlabel(r'${D} \log^2(D)$', fontsize=15)
    plt.ylabel('Sweep Times', fontsize=15)

    plt.legend(fontsize=14, loc='upper left', frameon=False)

    plt.tight_layout()

    plt.show()


def draw2():
    xvalue_20, yvalue_20 = read_width_file(
        '/data_toymodel_pulse/rule110/width_sweep', 20)
    xvalue_30, yvalue_30 = read_width_file(
        '/data_toymodel_pulse/rule110/width_sweep', 30)
    coefficients = np.polyfit(
        xvalue_20, yvalue_20, 1)  # 1 表示一次多项式（线性）
    slope1, intercept1 = coefficients
    y_fit_20 = slope1 * xvalue_20 + intercept1

    coefficients = np.polyfit(
        xvalue_30, yvalue_30, 1)  # 1 表示一次多项式（线性）
    slope2, intercept2 = coefficients
    y_fit_30 = slope2 * xvalue_30 + intercept2

    fig = plt.figure(dpi=200)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'Times New Roman'

    ax1 = plt.subplot(1, 1, 1)
    ax1.scatter(xvalue_20, yvalue_20, marker='o', edgecolors='#9a031e',
                facecolors='none', linewidth=3, s=60, label='Depth=20', zorder=5)
    ax1.plot(xvalue_20, y_fit_20, '--', c='#03045e',
             linewidth=2.5, zorder=5, label='Fit of Depth=20')
    ax1.scatter(xvalue_30, yvalue_30, marker='o', edgecolors='#fcbf49',
                facecolors='none', linewidth=3, s=60, label='Depth=30')
    ax1.plot(xvalue_30, y_fit_30, '-', c='#0096c7',
             linewidth=2.5, zorder=3, label='Fit of Depth=30')
    ax = plt.gca()
    mylinewidth = 1.5
    ax.spines['bottom'].set_linewidth(mylinewidth)
    ax.spines['left'].set_linewidth(mylinewidth)
    ax.spines['top'].set_linewidth(mylinewidth)
    ax.spines['right'].set_linewidth(mylinewidth)

    plt.xlim(left=5, right=50)
    plt.ylim(bottom=2000, top=10000)

    yticks_name = ['0'] + [f'{y}'+r'$\times 10^3$' for y in [
        2, 4, 6, 8]] + [f'{y}' + r'$\times 10^4$' for y in [1]]
    yticks = [y*1e3 for y in [0, 2, 4, 6, 8, 10]]
    # print(yticks_name)
    plt.xticks(fontsize=12)
    plt.yticks(yticks, yticks_name, fontsize=12, fontname='Times New Roman')

    plt.xlabel('Width', fontsize=15)
    plt.ylabel('Sweep Times', fontsize=15)

    plt.legend(fontsize=14, loc='lower right', frameon=False)

    plt.tight_layout()

    plt.show()


draw2()

# yticks = [0, 0.5, 1, 1.5, 2, 3]
# name_yticks = ['0', '1/2', '1', '3/2', '2', '2.5']
# xticks = [2*pi-2.5, 2*pi, 2*pi+2.5]
# # xticks = [1.4,1.5,1.6,1.7]

# # plt.yscale("log")
# plt.xlabel(r'$\Delta E_{\rm m} \tau$', fontsize=30)
# plt.ylabel(r'$\langle n \rangle_{\rm Con}$', fontsize=30)
# # plt.title('ibm_sherbrooke',fontsize = 40)
# # plt.yticks(yticks,name_yticks,fontsize=35)
# plt.yticks(yticks, name_yticks, fontsize=25)
# plt.xticks(xticks, (r'$2\pi-2.5$', r'$2\pi$', r'$2\pi+2.5$'), fontsize=25)
# # plt.xticks(xticks,('1.4','1.5','1.6','1.7'),fontsize=25)
# plt.xlim(2*pi-2.6, 2*pi + 2.6)
# plt.ylim(0, 2.5)


# plt.legend(fontsize=20, loc='lower right')

# ax2 = plt.subplot(1, 2, 2)

# ax2.plot(cp_tlist, cp_rm_theory(10, cp_tlist), '-',
#          c='#000000', linewidth=ll, zorder=3,  label='Theory')
# ax2.plot(cp_tlist, cp_rm_theory(20, cp_tlist),
#          '-', c='#000000', linewidth=ll, zorder=3)
# ax2.plot(cp_tlist, cp_rm_theory(40, cp_tlist),
#          '-', c='#000000', linewidth=ll, zorder=3)

# ax2.scatter(xrm10, yrm10, marker='o', edgecolors='#0077b6',
#             facecolors='none', linewidth=3, s=60, zorder=0, label=r"$L = 10$")
# ax2.scatter(xrm20, yrm20, marker='o', edgecolors='#2a9d8f',
#             facecolors='none', linewidth=3, s=60, label=r"$L = 20$")
# ax2.scatter(xrm40, yrm40, marker='o', edgecolors='#e76f51',
#             facecolors='none', linewidth=3, s=60, label=r"$L = 40$")

# yticks = [0, 0.5, 1, 1.5, 2, 3]
# name_yticks = ['0', '1/2', '1', '3/2', '2', '2.5']
# xticks = [2*pi-2.5, 2*pi, 2*pi+2.5]
# # xticks = [1.4,1.5,1.6,1.7]

# # plt.yscale("log")
# plt.xlabel(r'$\Delta E_{\rm m} \tau$', fontsize=30)
# plt.ylabel(r'$\langle n_R \rangle$', fontsize=30)
# # plt.title('ibm_sherbrooke',fontsize = 40)
# # plt.yticks(yticks,name_yticks,fontsize=35)
# plt.yticks(yticks, name_yticks, fontsize=25)
# plt.xticks(xticks, (r'$2\pi-2.5$', r'$2\pi$', r'$2\pi+2.5$'), fontsize=25)
# # plt.xticks(xticks,('1.4','1.5','1.6','1.7'),fontsize=25)
# plt.xlim(2*pi-2.6, 2*pi + 2.6)
# plt.ylim(0, 2.5)
# ax = plt.gca()
# ax.spines['bottom'].set_linewidth(2)
# ax.spines['left'].set_linewidth(2)
# ax.spines['top'].set_linewidth(2)
# ax.spines['right'].set_linewidth(2)

# plt.legend(fontsize=20, loc='lower right')

# plt.savefig('complete_graph.pdf')
# plt.show()
