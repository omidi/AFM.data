from imports import *
from toolbox import *


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


def align_two_cells(smaller, bigger):
    pad_width = 20
    mx = bigger[0]
    my = bigger[1]
    dx = smaller[0]
    dy = smaller[1]
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

def detect_troughs(cell, vicinity=None):
    if not vicinity:
        vicinity = int(np.ceil(max(cell[0])/.5))
    y = np.gradient(np.gradient(cell[1], 5), 5)  # 5 determines the smoothness of the profile
    mins = argrelextrema(y, np.greater, order=vicinity)
    troughs = []
    heights = []
    for m in mins[0]:
        if cell[0][m] < 1. or cell[0][m] > (max(cell[0] - 1.)):
            continue
        troughs.append(cell[0][m])
        heights.append((cell[1][m]))
    return troughs, heights


def plot_one_cell(cell):
    plt.plot(cell[0], cell[1], lw=2)
    p, h = detect_troughs(cell)
    plt.plot(p, h, 'o', color='black')
    plt.show()

def plot_derivativ_one_cell(cell):
    y = np.gradient(np.gradient(cell[1], 5), 5)
    plt.plot(cell[0], y, lw=2)
    p, h = detect_troughs(cell)
    plt.plot(p, [0 for i in xrange(len(p))], 'o', color='black')
    plt.show()


def find_daughter_cells(cell_name):
    cut_percentage = .05
    mother_last_time = load_file(return_all_files_for_cell(cell_name)[-1])
    daughter1 = map(load_file, return_all_files_for_cell(cell_name + '.1'))
    daughter2 = map(load_file, return_all_files_for_cell(cell_name + '.2'))
    plt.plot(mother_last_time[0], mother_last_time[1])
    troughs, heights = detect_troughs(mother_last_time)
    plt.plot(troughs, heights, 'o', color="red")
    o1, o2 = align_both_daughters(mother_last_time, daughter1[0], daughter2[0])
    t1, h1 = detect_troughs(daughter1[0])
    t2, h2 = detect_troughs(daughter2[0])
    plt.plot(daughter1[0][0] + o1, daughter1[0][1])
    plt.plot(daughter2[0][0] + o2, daughter2[0][1])
    plt.plot(t1 + o1, h1, 'o', color="black")
    plt.plot(t2 + o2, h2, 'o', color="black")
    plt.show()
    # offset1, offset2 = align_both_daughters(mother_last_time, daughter1[0], daughter2[0])
    # offset = align_two_cells(daughter1[1], daughter1[0]) + offset1
    # for i in xrange(1, 2):
    #     plt.plot(daughter1[i][0] + offset, daughter1[i][1])
    #     offset += align_two_cells(daughter1[i+1], daughter1[i]) + offset
    #
    # offset = align_two_cells(daughter2[1], daughter2[2]) + offset2
    # plt.plot(daughter2[1][0] + offset, daughter2[1][1])
    # plt.plot(daughter2[2][0], daughter2[2][1])
    # for i in xrange(1, 3):
    #     plt.plot(daughter2[i][0] + offset, daughter2[i][1])
    #     offset += align_two_cells(daughter2[i+1], daughter2[i]) + offset
    # plt.show()
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


def troughs_in_time(cellname):
    cells = map(load_file, return_all_files_for_cell(cellname))
    positions = [detect_troughs(cell)[0] for cell in cells]
    plt.xlim((-14, max(cells[-1][0]+14)))
    plt.ylim((-20, len(cells)+2))
    for t, c in enumerate(cells):
        x = len(cells) - t
        start = -(max(c[0]) - max(cells[-1][0]))/2
        plt.plot([start, start + max(c[0])], [x, x], lw=4)
    for t, pos in enumerate(positions):
        x = len(cells) - t
        start = -(max(cells[t][0]) - max(cells[-1][0]))/2
        plt.plot([p + start for p in pos], [x for i in xrange(len(pos))],
                 'o', color='black')
    d1 = map(load_file, return_all_files_for_cell(cellname + '.1'))
    d2 = map(load_file, return_all_files_for_cell(cellname + '.2'))
    o2 = align_two_cells(d2[0], cells[-1])
    o1 = align_two_cells(d1[0], cells[-1])
    for t in xrange(len(d1)):
        x = -t - 1
        start = o1 - (max(d1[t][0]) - max(d1[-1][0]))/2 - 5
        troughs = detect_troughs(d1[t])[0]
        plt.plot([start, start + max(d1[t][0])], [x, x], lw=4)
        plt.plot([t+start for t in troughs],
                 [x for i in xrange(len(troughs))], 'o', color='black')
    for t in xrange(len(d2)):
        x = -t - 1
        start = o2 - (max(d2[t][0]) - max(d2[-1][0]))/2 + 5
        plt.plot([start, start + max(d2[t][0])], [x, x], lw=4)
        troughs = detect_troughs(d2[t])[0]
        plt.plot([t+start for t in troughs],
                 [x for i in xrange(len(troughs))], 'o', color='black')
    plt.show()

# plot_all_cells_profile()
# find_daughter_cells('1.1')
# troughs_in_time('1.1')
plot_one_cell(load_file(return_all_files_for_cell('1.1')[1]))