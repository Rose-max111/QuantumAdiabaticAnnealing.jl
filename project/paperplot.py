from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import matplotlib.ticker as mticker


def read_depth_file(addpath, width):
    folder_path = os.getcwd() + addpath
    xvalue = []
    yvalue = []
    for filename in os.listdir(folder_path):
        # print(filename)
        if filename.startswith(f"W={width}_D="):
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
    # xvalue *= np.log(xvalue) * np.log(xvalue)
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

    xvalue_gradient, yvalue_gradient = read_depth_file('/data_toymodel/', 15)
    xvalue_pulse, yvalue_pulse = read_depth_file('/data_toymodel_pulse/', 15)

    # coefficients = np.polyfit(
    #     xvalue_gradient, yvalue_gradient, 1)  # 1 表示一次多项式（线性）
    # slope1, intercept1 = coefficients
    # y_fit_gradient = slope1 * xvalue_pulse + intercept1

    # coefficients = np.polyfit(xvalue_pulse, yvalue_pulse, 1)  # 1 表示一次多项式（线性）
    # slope2, intercept2 = coefficients
    # y_fit_pulse = slope2 * xvalue_pulse + intercept2

    # fig = plt.figure(dpi=200, figsize=(20, 5.5))
    fig = plt.figure(dpi=200)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'Times New Roman'

    ax1 = plt.subplot(1, 1, 1)
    ax1.scatter(xvalue_gradient, yvalue_gradient, marker='o', edgecolors='#9a031e',
                facecolors='none', linewidth=3, s=60, label='Energy Gradient', zorder=5)
    # ax1.plot(xvalue_pulse, y_fit_gradient, '--', c='#03045e',
    #          linewidth=2.5, zorder=5, label='Fit of Energy Gradient')
    ax1.scatter(xvalue_pulse, yvalue_pulse, marker='o', edgecolors='#fcbf49',
                facecolors='none', linewidth=3, s=60, label='Moving Heat Bath')
    # ax1.plot(xvalue_pulse, y_fit_pulse, '-', c='#0096c7',
    #          linewidth=2.5, zorder=3, label='Fit of Wave Packet')

    ax = plt.gca()
    mylinewidth = 1.5
    ax.spines['bottom'].set_linewidth(mylinewidth)
    ax.spines['left'].set_linewidth(mylinewidth)
    ax.spines['top'].set_linewidth(mylinewidth)
    ax.spines['right'].set_linewidth(mylinewidth)

    plt.xlim(left=0, right=34)
    plt.ylim(bottom=0, top=65000)

    yticks_name = ['0'] + [f'{y}'+r'$\times 10^4$' for y in [
        1, 2, 3, 4, 5, 6]]
    yticks = [y*1e4 for y in [0, 1, 2, 3, 4, 5, 6]]
    # print(yticks_name)
    plt.xticks(fontsize=15)
    plt.yticks(yticks, yticks_name, fontsize=15, fontname='Times New Roman')

    plt.xlabel(r'${M}$', fontsize=20)
    plt.ylabel(r'$T_c$', fontsize=20)

    plt.legend(fontsize=16, loc='upper left', frameon=False)

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
                facecolors='none', linewidth=3, s=60, label='T=20', zorder=5)
    ax1.plot(xvalue_20, y_fit_20, '--', c='#03045e',
             linewidth=2.5, zorder=5, label='Fit of T=20')
    ax1.scatter(xvalue_30, yvalue_30, marker='o', edgecolors='#fcbf49',
                facecolors='none', linewidth=3, s=60, label='T=30')
    ax1.plot(xvalue_30, y_fit_30, '-', c='#0096c7',
             linewidth=2.5, zorder=3, label='Fit of T=30')
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

    plt.xlabel('T', fontsize=15)
    plt.ylabel('Sweep Times', fontsize=15)

    plt.legend(fontsize=14, loc='lower right', frameon=False)

    plt.tight_layout()

    plt.show()


def modify(xdata, ydata):
    ydata = ydata / xdata
    xdata = np.log(xdata)
    return xdata, ydata


def draw3():
    fontsize_axis = 15
    fontsize_label = 20
    xdata1_5, ydata1_5 = read_depth_file('/data_toymodel_pulse/', 12)
    xdata1_5, ydata1_5 = modify(xdata1_5, ydata1_5)
    xdata1_3, ydata1_3 = read_depth_file('/data_toymodel_pulse/', 15)
    xdata1_3, ydata1_3 = modify(xdata1_3, ydata1_3)
    xdata1_24, ydata1_24 = read_depth_file('/data_toymodel_pulse/', 14)
    xdata1_24, ydata1_24 = modify(xdata1_24, ydata1_24)
    xdata1_18, ydata1_18 = read_depth_file('/data_toymodel_pulse/', 13)
    xdata1_18, ydata1_18 = modify(xdata1_18, ydata1_18)

    fig = plt.figure(dpi=200, figsize=(30, 5.5))
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'Times New Roman'

    ax1 = plt.subplot(1, 2, 1)
    ax1.scatter(xdata1_5, ydata1_5, marker='o', edgecolors='#273e47',
                facecolors='none', linewidth=3, s=60, label=r'$\Lambda = 1.5$', zorder=5)
    coefficients = np.polyfit(
        xdata1_5, ydata1_5, 1)  # 1 表示一次多项式（线性）
    slope1, intercept1 = coefficients
    ax1.scatter(xdata1_3, ydata1_3, marker='s', edgecolors='#bd632f',
                facecolors='none', linewidth=3, s=60, label=r'$\Lambda = 1.3$', zorder=5)
    coefficients = np.polyfit(
        xdata1_3, ydata1_3, 1)  # 1 表示一次多项式（线性）
    slope2, intercept2 = coefficients
    ax1.scatter(xdata1_24, ydata1_24, marker='^', edgecolors='#d8973c',
                facecolors='none', linewidth=3, s=60, label=r'$\Lambda = 1.24$', zorder=5)
    coefficients = np.polyfit(
        xdata1_24, ydata1_24, 1)  # 1 表示一次多项式（线性）
    slope3, intercept3 = coefficients
    ax1.scatter(xdata1_18, ydata1_18, marker='D', edgecolors='#d8c99b',
                facecolors='none', linewidth=3, s=60, label=r'$\Lambda = 1.18$', zorder=5)
    coefficients = np.polyfit(
        xdata1_18, ydata1_18, 1)  # 1 表示一次多项式（线性）
    slope4, intercept4 = coefficients

    # print(slope1, slope2, slope3, slope4)

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

    ax2 = plt.subplot(1, 2, 2)
    ax2.scatter(xvalue_20, yvalue_20, marker='o', edgecolors='#9a031e',
                facecolors='none', linewidth=3, s=60, label='M=20', zorder=5)
    ax2.plot(xvalue_20, y_fit_20, '--', c='#03045e',
             linewidth=2.5, zorder=5, label='Fit of M=20')
    ax2.scatter(xvalue_30, yvalue_30, marker='o', edgecolors='#fcbf49',
                facecolors='none', linewidth=3, s=60, label='M=30')
    ax2.plot(xvalue_30, y_fit_30, '-', c='#0096c7',
             linewidth=2.5, zorder=3, label='Fit of M=30')

    ax = plt.gca()
    mylinewidth = 1.5
    ax1.spines['bottom'].set_linewidth(mylinewidth)
    ax1.spines['left'].set_linewidth(mylinewidth)
    ax1.spines['top'].set_linewidth(mylinewidth)
    ax1.spines['right'].set_linewidth(mylinewidth)

    ax2.spines['bottom'].set_linewidth(mylinewidth)
    ax2.spines['left'].set_linewidth(mylinewidth)
    ax2.spines['top'].set_linewidth(mylinewidth)
    ax2.spines['right'].set_linewidth(mylinewidth)

    ax1.set_xlim(left=0.6, right=1.55)
    ax1.set_ylim(bottom=4, top=10)
    ax2.set_xlim(left=2, right=52)
    ax2.set_ylim(bottom=2000, top=10500)

    ax1.text(0.92, 0.98, '(a)', transform=ax1.transAxes,
             fontsize=18, verticalalignment='top')
    ax2.text(0.92, 0.98, '(b)', transform=ax2.transAxes,
             fontsize=18, verticalalignment='top')

    yticks_name = ['0'] + [f'{y}'+r'$\times 10^3$' for y in [
        2, 4, 6, 8]] + [f'{y}' + r'$\times 10^4$' for y in [1]]
    yticks = [y*1e3 for y in [0, 2, 4, 6, 8, 10]]
    # print(yticks_name)
    # ax2.set_xticks(fontsize=12)
    ax2.set_yticks(yticks, yticks_name, fontsize=fontsize_axis,
                   fontname='Times New Roman')
    ax2.tick_params(axis='x', labelsize=fontsize_axis)
    ax1.tick_params(axis='both', labelsize=fontsize_axis)

    # plt.ylim(bottom=0)
    # plt.xlim(left=0.6, right=1.55)
    # plt.ylim(bottom=4, top=10)
    ax1.set_xlabel(r'$\ln{(\ln(M))}$', fontsize=fontsize_label)
    ax1.set_ylabel(
        r'$\ln\left(v\right)$', fontsize=fontsize_label)

    ax2.set_xlabel(r'$N$', fontsize=fontsize_label)
    ax2.set_ylabel(r'$T_c$', fontsize=fontsize_label)
    ax1.legend(fontsize=16, loc='upper left', frameon=False)
    ax2.legend(fontsize=16, loc='upper left', frameon=False)

    # plt.tight_layout()
    plt.subplots_adjust(left=0.08, right=0.98, top=0.98,
                        bottom=0.11, wspace=0.23, hspace=0.2)

    plt.savefig('main_result.pdf', format='pdf')
    plt.show()


def draw4():
    fontsize_axis = 20
    fontsize_label = 25
    xdata1_5, ydata1_5 = read_depth_file('/data_toymodel_pulse/', 12)
    xdata1_5, ydata1_5 = modify(xdata1_5, ydata1_5)
    xdata1_3, ydata1_3 = read_depth_file('/data_toymodel_pulse/', 15)
    xdata1_3, ydata1_3 = modify(xdata1_3, ydata1_3)
    xdata1_24, ydata1_24 = read_depth_file('/data_toymodel_pulse/', 14)
    xdata1_24, ydata1_24 = modify(xdata1_24, ydata1_24)
    xdata1_18, ydata1_18 = read_depth_file('/data_toymodel_pulse/', 13)
    xdata1_18, ydata1_18 = modify(xdata1_18, ydata1_18)

    fig = plt.figure(dpi=200, figsize=(20, 6.5))
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'Times New Roman'

    ax1 = plt.subplot(1, 3, 1)
    ax1.scatter(xdata1_5, ydata1_5, marker='o', edgecolors='#273e47',
                facecolors='none', linewidth=3, s=60, label=r'$\Lambda = 1.5$', zorder=5)
    coefficients = np.polyfit(
        xdata1_5, ydata1_5, 1)  # 1 表示一次多项式（线性）
    slope1, intercept1 = coefficients
    ax1.scatter(xdata1_3, ydata1_3, marker='s', edgecolors='#bd632f',
                facecolors='none', linewidth=3, s=60, label=r'$\Lambda = 1.3$', zorder=5)
    coefficients = np.polyfit(
        xdata1_3, ydata1_3, 1)  # 1 表示一次多项式（线性）
    slope2, intercept2 = coefficients
    ax1.scatter(xdata1_24, ydata1_24, marker='^', edgecolors='#d8973c',
                facecolors='none', linewidth=3, s=60, label=r'$\Lambda = 1.24$', zorder=5)
    coefficients = np.polyfit(
        xdata1_24, ydata1_24, 1)  # 1 表示一次多项式（线性）
    slope3, intercept3 = coefficients
    ax1.scatter(xdata1_18, ydata1_18, marker='D', edgecolors='#d8c99b',
                facecolors='none', linewidth=3, s=60, label=r'$\Lambda = 1.18$', zorder=5)
    coefficients = np.polyfit(
        xdata1_18, ydata1_18, 1)  # 1 表示一次多项式（线性）
    slope4, intercept4 = coefficients

    # print(slope1, slope2, slope3, slope4)

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

    ax2 = plt.subplot(1, 3, 2)
    ax2.scatter(xvalue_20, yvalue_20, marker='o', edgecolors='#9a031e',
                facecolors='none', linewidth=3, s=60, label='M=20', zorder=5)
    ax2.plot(xvalue_20, y_fit_20, '--', c='#03045e',
             linewidth=2.5, zorder=5, label='Fit of M=20')
    ax2.scatter(xvalue_30, yvalue_30, marker='o', edgecolors='#fcbf49',
                facecolors='none', linewidth=3, s=60, label='M=30')
    ax2.plot(xvalue_30, y_fit_30, '-', c='#0096c7',
             linewidth=2.5, zorder=3, label='Fit of M=30')

    ax = plt.gca()
    mylinewidth = 1.5
    ax1.spines['bottom'].set_linewidth(mylinewidth)
    ax1.spines['left'].set_linewidth(mylinewidth)
    ax1.spines['top'].set_linewidth(mylinewidth)
    ax1.spines['right'].set_linewidth(mylinewidth)

    ax2.spines['bottom'].set_linewidth(mylinewidth)
    ax2.spines['left'].set_linewidth(mylinewidth)
    ax2.spines['top'].set_linewidth(mylinewidth)
    ax2.spines['right'].set_linewidth(mylinewidth)

    # ax1.set_ylim(bottom=4, top=10)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlim(left=1.85)
    ax2.set_xlim(left=2, right=52)
    ax2.set_ylim(bottom=2000, top=10500)

    ax1.text(0.92, 0.98, '(a)', transform=ax1.transAxes,
             fontsize=18, verticalalignment='top')
    ax2.text(0.92, 0.98, '(b)', transform=ax2.transAxes,
             fontsize=18, verticalalignment='top')

    yticks_name = ['0'] + [f'{y}'+r'$\times 10^3$' for y in [
        2, 4, 6, 8]] + [f'{y}' + r'$\times 10^4$' for y in [1]]
    yticks = [y*1e3 for y in [0, 2, 4, 6, 8, 10]]
    # print(yticks_name)
    # ax2.set_xticks(fontsize=12)
    ax2.set_yticks(yticks, yticks_name, fontsize=fontsize_axis,
                   fontname='Times New Roman')
    ax2.tick_params(axis='x', labelsize=fontsize_axis)
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_axis)

    ax1.tick_params(axis='x', which='minor', labelsize=fontsize_axis)

    # plt.ylim(bottom=0)
    # plt.xlim(left=0.6, right=1.55)
    # plt.ylim(bottom=4, top=10)
    ax1.set_xlabel(r'$\ln(M)$', fontsize=fontsize_label)
    ax1.set_ylabel(
        r'$v$', fontsize=fontsize_label)

    ax2.set_xlabel(r'$N$', fontsize=fontsize_label)
    ax2.set_ylabel(r'$T_c$', fontsize=fontsize_label)
    ax1.legend(fontsize=16, loc='upper left', frameon=False)
    ax2.legend(fontsize=16, loc='upper left', frameon=False)

    # plt.tight_layout()
    # plt.subplots_adjust(left=0.08, right=0.98, top=0.98,
    #                     bottom=0.11, wspace=0.23, hspace=0.2)

    xvalue_pulse, yvalue_pulse = read_depth_file('/data_toymodel_pulse/', 15)
    xvalue_fixinput, yvalue_fixinput = read_depth_file(
        '/data_toymodel_fixinput/', 15)
    xvalue_fixoutput, yvalue_fixoutput = read_depth_file(
        '/data_toymodel_fixoutput/', 15)
    xvalue_reverse, yvalue_reverse = read_depth_file(
        '/data_toymodel_reverse_pulse/', 12)
    # yvalue_pulse = np.log(yvalue_pulse)
    # yvalue_fixinput = np.log(yvalue_fixinput)
    # yvalue_fixoutput = np.log(yvalue_fixoutput)
    # yvalue_reverse = np.log(yvalue_reverse)

    ax3 = plt.subplot(1, 3, 3)
    ax3.scatter(xvalue_pulse, yvalue_pulse, marker='o', edgecolors='#273e47',
                facecolors='none', linewidth=3, s=60, label='HB(Forward)', zorder=5)
    ax3.scatter(xvalue_fixinput, yvalue_fixinput, marker='s', edgecolors='#bd632f',
                facecolors='none', linewidth=3, s=60, label='GA(Fixed input)', zorder=5)
    ax3.scatter(xvalue_reverse, yvalue_reverse, marker='^', edgecolors='#d8973c',
                facecolors='none', linewidth=3, s=60, label='HB(Backward)', zorder=5)
    ax3.scatter(xvalue_fixoutput, yvalue_fixoutput, marker='D', edgecolors='#d8c99b',
                facecolors='none', linewidth=3, s=60, label='GA(Fixed output)', zorder=5)
    ax3.text(0.92, 0.98, '(c)', transform=ax3.transAxes,
             fontsize=18, verticalalignment='top')

    xvalue_line = np.linspace(1.0, 16, 100)
    yvalue_linear = np.exp(xvalue_line) * 280
    yvalue_log = xvalue_line * \
        np.power(np.log(xvalue_line), 2.27) * 70

    ax3.plot(xvalue_line, yvalue_linear, '--', c='#03045e',
             linewidth=2.5, zorder=5, label=r'$y \propto e^x$')
    ax3.plot(xvalue_line, yvalue_log, '-', c='#0096c7',
             linewidth=2.5, zorder=3, label=r'$y \propto x\ln^{2.27}(x)$')

    mylinewidth = 1.5
    ax3.spines['bottom'].set_linewidth(mylinewidth)
    ax3.spines['left'].set_linewidth(mylinewidth)
    ax3.spines['top'].set_linewidth(mylinewidth)
    ax3.spines['right'].set_linewidth(mylinewidth)
    ax3.set_xlim(left=1.0, right=13)
    # ax3.set_ylim(bottom=1.7, top=20)
    ax3.set_yscale('log')
    ax3.set_ylim(bottom=1e2, top=3e7)

    ax3.tick_params(axis='both', labelsize=fontsize_axis)

    ax3.set_xlabel(r'$M$', fontsize=fontsize_label)
    ax3.set_ylabel(
        r'$T_c$', fontsize=fontsize_label)
    ax3.legend(fontsize=16, loc='upper left', frameon=False)

    plt.subplots_adjust(left=0.08, right=0.98, top=0.98,
                        bottom=0.11, wspace=0.3, hspace=0.2)

    plt.savefig(
        '/Users/rose/gitrepo/TCnotes/notes/images/main_result.pdf', format='pdf')
    # plt.show()


# draw3()
# draw1()
# draw2()
draw4()

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
