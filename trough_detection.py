from imports import *
from toolbox import *


def detect_troughs(signal, vicinity=10):
    y = np.gradient(np.gradient(signal, 5), 5)  # 5 determines the smoothness of the profile
    mins = argrelextrema(y, np.greater, order=vicinity)
    return mins[0]


def plot_one_cell_in_time(cell_name):
    data = map(load_file, return_all_files_for_cell(cell_name))
    orders_by_length = np.argsort(([max(d[0]) for d in data]))[::-1]
    longest_cell = orders_by_length[0]
    plt.figure(figsize=(12,12))
    plt.plot(data[longest_cell][0], data[longest_cell][1] / max(data[longest_cell][1]))
    plt.ylim((0,1.2))
    total_pad = 0.
    for t in xrange(1, len(data)):
        t_prev = orders_by_length[t-1]
        t_curr = orders_by_length[t]
        x = data[t_curr][1]
        y = data[t_prev][1]
        trough_indices = detect_troughs(x)
        conv = np.convolve(y, x[::-1])
        left_pad = np.argmax(conv) - len(x) + 1
        if left_pad<0:
            left_pad = 0
        total_pad += data[t_prev][0][left_pad]
        plt.plot(data[t_curr][0] + total_pad, (x / max(x)))
        plt.plot(np.take(data[t_curr][0], trough_indices) + total_pad,
                 np.take(x, trough_indices) / max(x), 'D', color='black')
    plt.show()


def plot_all_cells_profile():
    all_files = return_all_files()
    plt.ylim((0.2,1.1))
    plt.xlim((-.1,1.1))
    plt.xlabel("Normalized cell length")
    plt.ylabel("Normalized surface height")
    x = []
    y = []
    for a_file in all_files:
        data = load_file(a_file)
        plt.plot(data[0]/max(data[0]), data[1]/max(data[1]), '-', color='black', alpha=.01)
        x.append(data[0]/max(data[0]))
        y.append(data[1]/max(data[1]))
    x_sm = np.array([item for sublist in x for item in sublist])
    y_sm = np.array([item for sublist in y for item in sublist])
    f = UnivariateSpline(x_sm, y_sm, k=5)
    x_smooth = np.linspace(0., 1., 100)
    y_smooth = f(x_smooth)
    plt.plot(x_smooth, y_smooth, color='red', lw=3)
    plt.show()


def align_two_cells(mother, daughter):
    """
    it is assumed that the mother cell is longer than the daughter cell
    :param mother:
    :param daughter:
    :return:
    """
    pad_width = 20
    mx = mother[0]
    my = mother[1]
    dx = daughter[0]
    dy = daughter[1]
    dyd = np.gradient(dy)  # 5 determines the smoothness of the profile
    myd = np.gradient(my)  # 5 determines the smoothness of the profile
    myd = np.pad(myd, (pad_width, pad_width), 'constant')
    distance = []
    distance_rev = []
    for pad in xrange(len(myd) - len(dyd) + 1):
        dist = np.mean(([(p[0] - p[1])**2 for p in zip(dyd, myd[pad:(pad+len(dyd))])]))
        dist_rev = np.mean(([(p[0] - p[1])**2 for p in zip(dyd[::-1], myd[pad:(pad+len(dyd))])]))
        distance.append(dist)
        distance_rev.append(dist_rev)
    if min(distance) < min(distance_rev):
        left_pad = np.argmin(distance) - pad_width
        if left_pad > 0:
            offset = mx[left_pad]
        else:
            offset = -mx[-left_pad]
    else:
        left_pad = np.argmin(distance_rev) - pad_width
        if left_pad > 0:
            offset = mx[left_pad]
        else:
            offset = -mx[-left_pad]
    return offset


def align_both_daughters(mother, daughter1, daughter2):
    plt.plot(mother[0], mother[1])
    offset1 = align_two_cells(mother, daughter1)
    plt.plot(daughter1[0] + offset1, daughter1[1])
    offset2 = align_two_cells(mother, daughter2)
    plt.plot(daughter2[0] + offset2, daughter2[1])
    return offset1, offset2


def find_daughter_cells(cell_name):
    cut_percentage = .05
    mother_last_time = load_file(return_all_files_for_cell(cell_name)[-1])
    daughter1 = map(load_file, return_all_files_for_cell(cell_name + '.1'))
    daughter2 = map(load_file, return_all_files_for_cell(cell_name + '.2'))
    offset1, offset2 = align_both_daughters(mother_last_time, daughter1[0], daughter2[0])
    # offset = align_two_cells(daughter1[1], daughter1[0]) + offset1
    # for i in xrange(1, 2):
    #     plt.plot(daughter1[i][0] + offset, daughter1[i][1])
    #     offset += align_two_cells(daughter1[i+1], daughter1[i]) + offset
    #
    offset = align_two_cells(daughter2[1], daughter2[2]) + offset2
    plt.plot(daughter2[1][0] + offset, daughter2[1][1])
    plt.plot(daughter2[2][0], daughter2[2][1])
    for i in xrange(1, 3):
        plt.plot(daughter2[i][0] + offset, daughter2[i][1])
        offset += align_two_cells(daughter2[i+1], daughter2[i]) + offset
    plt.show()
    exit()
    cut = int(len(mother_last_time[0])*cut_percentage)
    mother_x = mother_last_time[0][cut:(len(mother_last_time[0]) - cut)]
    mother_y = mother_last_time[1][cut:(len(mother_last_time[1]) - cut)]
    plt.plot(mother_x, mother_y / max(mother_y))
    total_pad = 0
    y = mother_y
    x = daughter1[0][1]
    cut = int(len(x)*cut_percentage)
    x = x[cut:(len(x)-cut)]
    xx = np.gradient(x, 5)  # 5 determines the smoothness of the profile
    yy = np.gradient(y, 5)  # 5 determines the smoothness of the profile
    conv = np.convolve(yy, xx[::-1])
    left_pad = np.argmax(conv) - len(x) + 1
    if left_pad<0:
        left_pad = -1 * daughter1[0][0][abs(left_pad)]
    total_pad += mother_last_time[0][left_pad]
    x_axis = daughter1[0][0]
    x_axis = x_axis[cut:(len(x_axis)-cut)]
    plt.plot(x_axis + total_pad, (x / max(x)))
    total_pad = 0
    y = mother_y
    x = daughter2[0][1]
    cut = int(len(x)*cut_percentage)
    x = x[cut:(len(x)-cut)]
    xx = np.gradient(x, 5)  # 5 determines the smoothness of the profile
    yy = np.gradient(y, 5)  # 5 determines the smoothness of the profile
    conv = np.convolve(yy, xx[::-1])
    left_pad = np.argmax(conv) - len(x) + 1
    if left_pad<0:
        left_pad = -1 * daughter2[0][0][abs(left_pad)]
    total_pad += mother_last_time[0][left_pad]
    x_axis = daughter2[0][0]
    x_axis = x_axis[cut:(len(x_axis)-cut)]
    plt.plot(x_axis + total_pad, (x / max(x)))
    plt.show()


# plot_all_cells_profile()
find_daughter_cells('1.1.1')
