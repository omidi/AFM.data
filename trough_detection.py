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



def align_daughter(mother, daughter):
    """
    It is assumed that the length of daughter is shorter than mother
    :param mother: mother data as given by load_file() function
    :param daughter: daughter data
    :return: is the offset of daughter and x-axis position relative to its mother
    """
    cut_percentage = .05
    cut = int(len(mother[0])*cut_percentage)
    mx = mother[0][cut:(len(mother[0]) - cut)]
    my = mother[1][cut:(len(mother[1]) - cut)]

    cut = int(len(daughter[0])*cut_percentage)
    dx = daughter[0][cut:(len(daughter[0]) - cut)]
    dy = daughter[1][cut:(len(daughter[1]) - cut)]

    dyd = np.gradient(dy, 5)  # 5 determines the smoothness of the profile
    myd = np.gradient(my, 5)  # 5 determines the smoothness of the profile

    plt.plot(mx, my)
    # plt.plot(dx, dy)
    # plt.show()
    # exit()
    best_score = np.inf
    best_offset = 0
    final_dy = np.array([])
    for rev in [False, True]:
        if rev:
            y = dy[::-1]
            h = dy[::-1]
        else:
            y = dy
            h = dy
        conv = np.convolve(my, y)
        # plt.plot(conv)
        left_pad = np.argmax(conv) - len(y) + 1
        # print left_pad, np.argmax(conv)
        # continue
        if left_pad<0:
            offset = -mx[abs(left_pad)] if not rev else mx[-abs(left_pad)] - max(mx)
            y_padded = y
        else:
            offset = mx[left_pad] if not rev else max(mx) - mx[-left_pad]
            y_padded = np.pad(y, (left_pad, 0), 'constant')
        # plt.plot(dx + offset, h)
        # plt.show()
        # exit()

        # dx_padded = np.pad(dx,(left_pad, right_pad), 'edge')
        score = np.sum(np.power([np.diff(p) for p in zip(y_padded, my)], 2))
        if score < best_score:
            best_offset = offset
            final_y = h
            best_score = score
    plt.plot(dx + best_offset, final_y)
    plt.show()


def align_daughter2(mother, daughter):
    cut_percentage = .0
    cut = int(len(mother[0])*cut_percentage)
    mx = mother[0][cut:(len(mother[0]) - cut)]
    my = mother[1][cut:(len(mother[1]) - cut)]
    cut = int(len(daughter[0])*cut_percentage)
    dx = daughter[0][cut:(len(daughter[0]) - cut)]
    dy = daughter[1][cut:(len(daughter[1]) - cut)]
    pad_width = len(dy)/2
    dyd = np.gradient(dy, 5)  # 5 determines the smoothness of the profile
    myd = np.gradient(my, 5)  # 5 determines the smoothness of the profile
    my_pad = np.pad(myd, (len(dyd)-pad_width, len(dyd)+pad_width), 'constant')

    min_dist = np.inf
    best_padding = 0
    total_padding = len(my_pad) - len(dy) + 1
    inverse = False
    for padding in xrange(total_padding):
        dist = np.sum(([(p[0] - p[1])**2 for p in zip(dyd, my_pad[padding:(padding+len(dy))])]))
        if dist < min_dist:
            min_dist = dist
            best_padding = padding

    print best_padding, total_padding, min_dist, pad_width
    # for padding in xrange(total_padding):
    #     dist = np.sum(([(p[0] - p[1])**2 for p in zip(dy[::-1], my_pad[padding:(padding+len(dy))])]))
        # if dist < min_dist:
        #     min_dist = dist
        #     best_padding = padding
        #     inverse = True

    if best_padding < pad_width:
        offset = -dx[best_padding]
    else:
        offset = mx[best_padding - pad_width]

    plt.plot(mx, my)
    if not inverse:
        plt.plot(dx + offset, dy)
    else:
        plt.plot(dx + offset, dy[::-1])
    plt.show()
    print best_padding, total_padding, offset, min_dist, inverse


def find_daughter_cells(cell_name):
    cut_percentage = .05
    mother_last_time = load_file(return_all_files_for_cell(cell_name)[-1])
    daughter1 = map(load_file, return_all_files_for_cell(cell_name + '.1'))
    daughter2 = map(load_file, return_all_files_for_cell(cell_name + '.2'))
    align_daughter2(mother_last_time, daughter2[0])
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
