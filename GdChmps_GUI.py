import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as transforms
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as sg


def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


def overlap(imag, iv, rotd):
    magx = 114/imag
    magy = 85.7/imag
    su = magx * magy

    iv_rot = [iv[0] * np.cos(rotd * np.pi / 180) - iv[1] * np.sin(rotd * np.pi / 180),
              iv[0] * np.sin(rotd * np.pi / 180) + iv[1] * np.cos(rotd * np.pi / 180)]

    if (np.abs(iv_rot[0]) > magx) | (np.abs(iv_rot[1]) > magy):
        si = 0
    else:
        si = np.min([np.abs((iv_rot[0] - magx) * (iv_rot[1] + magy)),
                     np.abs((iv_rot[0] + magx) * (iv_rot[1] - magy)),
                     np.abs((iv_rot[0] - magx) * (iv_rot[1] - magy)),
                     np.abs((iv_rot[0] + magx) * (iv_rot[1] + magy))])

    ol = 100 * si / su  # overlap in %
    return ol


def draw_rectangle(vxo, vyo, ind1, ind2, rotd, imag, ax):
    t_r = transforms.Affine2D().rotate_deg(rotd)
    for k in range(ind1, ind2):
        rect = patches.Rectangle((0, 0), 114 / imag, 85.7 / imag,
                                 linewidth=1,
                                 edgecolor='k',
                                 facecolor='none',
                                 alpha=0.5)

        vcx = 114 / imag / 2
        vcy = 85.7 / imag / 2
        vcx_p = vcx * np.cos(rotd * np.pi / 180) - vcy * np.sin(rotd * np.pi / 180)
        vcy_p = vcx * np.sin(rotd * np.pi / 180) + vcy * np.cos(rotd * np.pi / 180)
        t_t = transforms.Affine2D().translate(1000 * vxo[k] - vcx_p,
                                              1000 * vyo[k] - vcy_p)
        rect.set_transform(t_r + t_t + ax.transData)

        ax.add_patch(rect)


def canvas_layout(nc):
    clayout = [[sg.Canvas(key='olCanvas'+str(nc))]]
    return clayout


def open_ol_window(ntabs):
    layout2 = [
        [sg.TabGroup([[sg.Tab('S' + str(i + 1).zfill(3), canvas_layout(i)) for i in range(0, ntabs)]]),
         sg.VSeparator(),
         sg.Column([[sg.Text('Overlap :\n==============')],
                    [sg.Text('', key='output')]])]
    ]

    windw = sg.Window('Overlap information', layout2, finalize=True, resizable=True,
                      element_justification="left")
    return windw


# column1 = [
#     [sg.Text('Set beam shift to 0 before taking initial and end points.')],
#     [sg.Text('Tilt and Rotation must not change between the initial and final points of a serie.')],
#     [sg.Text('Initial points:', tooltip='Point Lists are saved in C:\\ProgramData\\Carl Zeiss\\SmartSEM\\User\\sem\n'
#                                         'A line scan must be defined by two consecutive points\n'
#                                         'A map must be defined by 4 consecutive points',
#              size=(20, 1))],
#     [sg.Input(key='ptInit', size=(60, 1)), sg.FileBrowse(target='ptInit', size=(6, 1),
#                                                          initial_folder='C:/ProgramData/Carl Zeiss/SmartSEM/User/sem')],
#     [sg.Text('# images of each series: ',
#              tooltip='Enter the number of images of each series as an integer separated by a comma\n'
#                      'e.g. 2,3,4 for 3 series of 2, 3, and 4 images respectively.\n'
#                      'For maps, enter m, the # of images in the x-axis, and n, the # of images in the y-axis, as mxn',
#              size=(20, 1)), sg.Input(key='n_im')],
#     [sg.Text('Magnification of each series: ', tooltip='Enter the magnification of each series separated by a comma\n'
#                                                        'e.g. 1000,3000,2000 for 3 series with mag = 1k, 3k, and 2k '
#                                                        'respectively.',
#              size=(20, 1)), sg.Input(key='mag')],
#     [sg.Text('Scan rot :', tooltip='Enter the scan rotation value',
#              size=(20, 1)), sg.Input(key='rot', default_text='0')],
#     [sg.Checkbox('Display overlap information', default=False, key='ol_cb')],
#     [sg.Text('Save as: ', size=(20, 1))],
#     [sg.Input(key='saveName', size=(58, 1)), sg.SaveAs(target='saveName', size=(8, 1))],
#     [sg.Text('In Stage Points List > On Goto, check Move XY Only.')],
#     [sg.Submit('Generate point list'), sg.Exit()]
# ]
column1 = [
    [sg.Text('Set beam shift to 0 before taking initial and end points.')],
    [sg.Text('Tilt and Rotation must not change between the initial and final points of a serie.')],
    [sg.Text('Initial points:', tooltip='Point Lists are saved in C:\\ProgramData\\Carl Zeiss\\SmartSEM\\User\\sem\n'
                                        'A line scan must be defined by two consecutive points\n'
                                        'A map must be defined by 4 consecutive points',
             size=(20, 1))],
    [sg.Input(key='ptInit', size=(60, 1)),
     sg.FileBrowse(target='ptInit',
                   size=(6, 1),
                   initial_folder='D:/_Ressources/SEM/GdChamps_v1.x/py/IN_2_4_2.TXT')],
    [sg.Text('# images of each series: ',
             tooltip='Enter the number of images of each series as an integer separated by a comma\n'
                     'e.g. 2,3,4 for 3 series of 2, 3, and 4 images respectively.\n'
                     'For maps, enter m, the # of images in the x-axis, and n, the # of images in the y-axis, as mxn',
             size=(20, 1)), sg.Input(key='n_im',
                                     default_text='5,5x4,8')],
    [sg.Text('Magnification of each series: ', tooltip='Enter the magnification of each series separated by a comma\n'
                                                       'e.g. 1000,3000,2000 for 3 series with mag = 1k, 3k, and 2k '
                                                       'respectively.',
             size=(20, 1)), sg.Input(key='mag', default_text='10,10,10')],
    [sg.Text('Scan rot :', tooltip='Enter the scan rotation value',
             size=(20, 1)), sg.Input(key='rot', default_text='60')],
    [sg.Checkbox('Display overlap information', default=False, key='ol_cb')],
    [sg.Text('Save as: ', size=(20, 1))],
    [sg.Input(key='saveName', size=(58, 1)), sg.SaveAs(target='saveName', size=(8, 1))],
    [sg.Text('In Stage Points List > On Goto, check Move XY Only.')],
    [sg.Submit('Generate point list'), sg.Exit()]
]
column2 = [
    [sg.Canvas(key='figCanvas')],
    # [sg.Text('Overlap :\n==============')],
    # [sg.Text('', key='output')],
    [sg.Text('Estimated time :\n==============')],
    [sg.Text('',
             size=(40, 1), key='time')]
]
layout = [
    [sg.Column(column1),
     sg.VSeparator(),
     sg.Column(column2)]
]

fig = plt.figure()
ax1 = fig.gca()

ax1.plot(65, 65, 'kx')
c1 = plt.Circle((65, 65), 25, color='k', fill=False)
c2 = plt.Circle((65, 65), 40, color='k', fill=False)
ax1.add_artist(c1)
ax1.add_artist(c2)

ax1.set_xlim([25, 105])
ax1.set_ylim([25, 105])
ax1.set_aspect('equal')
ax1.set_xlabel('X [mm]')
ax1.set_ylabel('Y [mm]')

window = sg.Window('Point List Generator', layout, finalize=True, resizable=True, element_justification="left")
fig_c_a = draw_figure(window['figCanvas'].TKCanvas, fig)

while True:
    event, values = window.read()
    if event in (None, 'Exit'):
        break

    ptListIn = values['ptInit']
    ptListOut = values['saveName']

    n_im = []
    mag = []
    rot = 0
    while True:
        fig_c_a.get_tk_widget().forget()
        plt.clf()

        try:
            n_im = values['n_im'].split(',')
            n_im = [i.split('x') for i in n_im]
            n_im = [[int(i) for i in j] for j in n_im]
            mag = [float(i) for i in values['mag'].split(',')]
            rot = values['rot']
        except ValueError:
            sg.Popup('Error', 'Please enter # of images as integers (or in the form of mxn) separated by a comma\n'
                              'Please enter magnification as doubles separated by a comma')
            event, values = window.read()
            continue
        break

    with open(ptListIn, 'r') as pl:
        pll = pl.readlines()

        v_x = []
        v_y = []
        v_z = []
        v_t = []
        v_r = []
        v_wd = []
        for i in range(4, len(pll)):
            tmp = pll[i].split(',')
            v_x.append(float(tmp[1]))
            v_y.append(float(tmp[2]))
            v_z.append(float(tmp[3]))
            v_t.append(float(tmp[4]))
            v_r.append(float(tmp[5]))
            v_wd.append(float(tmp[8][0:-1]))

        lnim = len([i for j in n_im for i in j])
        totim = np.sum([np.prod(j) for j in n_im])

    if (v_t[0] == v_t[-1]) & (v_r[0] == v_r[-1]):
        if len(v_x) == 2*lnim:
            fig = plt.figure()
            ax1 = fig.gca()

            ax1.plot(65, 65, 'kx')
            c1 = plt.Circle((65, 65), 25, color='k', fill=False)
            c2 = plt.Circle((65, 65), 40, color='k', fill=False)
            ax1.add_artist(c1)
            ax1.add_artist(c2)

            with open(ptListOut, 'w') as pl:
                pl.write('Leo Points List\n')
                pl.write('Absolute\n')
                pl.write('Label,X,Y,Z,T,R,M,Mag,WD\n')
                pl.write(str(totim + 1) + '\n')

                c = 0
                v_x_o = []
                v_y_o = []
                v_z_o = []
                v_wd_o = []

                l_ol = []
                ol_fig = []
                ol_ax = []
                for i in range(0, len(n_im)):
                    if len(n_im[i]) == 1:
                        # line A to B
                        lv_x_o1 = len(v_x_o)

                        v_x_o = np.concatenate((v_x_o, np.linspace(v_x[c], v_x[c + 1], n_im[i][0])), axis=0)
                        v_y_o = np.concatenate((v_y_o, np.linspace(v_y[c], v_y[c + 1], n_im[i][0])), axis=0)
                        v_z_o = np.concatenate((v_z_o, np.linspace(v_z[c], v_z[c + 1], n_im[i][0])), axis=0)
                        v_wd_o = np.concatenate((v_wd_o, np.linspace(v_wd[c], v_wd[c + 1], n_im[i][0])), axis=0)
                        c = c + 2

                        lv_x_o2 = len(v_x_o)
                        ax1.plot(1000 * v_x_o[lv_x_o1:lv_x_o2], 1000 * v_y_o[lv_x_o1:lv_x_o2], 'x-')

                        l_ol.append(overlap(mag[i],
                                            [1000 * (v_x_o[lv_x_o1 + 1] - v_x_o[lv_x_o1]),
                                             1000 * (v_y_o[lv_x_o1 + 1] - v_y_o[lv_x_o1])],
                                            int(rot)))

                    elif len(n_im[i]) == 2:
                        # trapezoid delimited by A B C D
                        lv_x_o1 = len(v_x_o)

                        AD_x = np.linspace(v_x[c], v_x[c + 3], n_im[i][1])
                        BC_x = np.linspace(v_x[c + 1], v_x[c + 2], n_im[i][1])
                        AD_y = np.linspace(v_y[c], v_y[c + 3], n_im[i][1])
                        BC_y = np.linspace(v_y[c + 1], v_y[c + 2], n_im[i][1])
                        AD_z = np.linspace(v_z[c], v_z[c + 3], n_im[i][1])
                        BC_z = np.linspace(v_z[c + 1], v_z[c + 2], n_im[i][1])
                        AD_wd = np.linspace(v_wd[c], v_wd[c + 3], n_im[i][1])
                        BC_wd = np.linspace(v_wd[c + 1], v_wd[c + 2], n_im[i][1])

                        for j in range(0, n_im[i][1]):
                            if (j % 2) == 1:
                                v_x_o = np.concatenate((v_x_o, np.linspace(AD_x[j], BC_x[j], n_im[i][0])), axis=0)
                                v_y_o = np.concatenate((v_y_o, np.linspace(AD_y[j], BC_y[j], n_im[i][0])), axis=0)
                                v_z_o = np.concatenate((v_z_o, np.linspace(AD_z[j], BC_z[j], n_im[i][0])), axis=0)
                                v_wd_o = np.concatenate((v_wd_o, np.linspace(AD_wd[j], BC_wd[j], n_im[i][0])), axis=0)
                            else:
                                v_x_o = np.concatenate((v_x_o, np.linspace(BC_x[j], AD_x[j], n_im[i][0])), axis=0)
                                v_y_o = np.concatenate((v_y_o, np.linspace(BC_y[j], AD_y[j], n_im[i][0])), axis=0)
                                v_z_o = np.concatenate((v_z_o, np.linspace(BC_z[j], AD_z[j], n_im[i][0])), axis=0)
                                v_wd_o = np.concatenate((v_wd_o, np.linspace(BC_wd[j], AD_wd[j], n_im[i][0])), axis=0)
                        c = c + 4

                        lv_x_o2 = len(v_x_o)
                        ax1.plot(1000 * v_x_o[lv_x_o1:lv_x_o2], 1000 * v_y_o[lv_x_o1:lv_x_o2], 'x-')

                        l_olx = overlap(mag[i],
                                        [1000 * (v_x_o[lv_x_o1] - v_x_o[lv_x_o1 + 2*n_im[i][0]-1]),
                                         1000 * (v_y_o[lv_x_o1] - v_y_o[lv_x_o1 + 2*n_im[i][0]-1])],
                                        int(rot))
                        l_oly = overlap(mag[i],
                                        [1000 * (v_x_o[lv_x_o1] - v_x_o[lv_x_o1 + 1]),
                                         1000 * (v_y_o[lv_x_o1] - v_y_o[lv_x_o1 + 1])],
                                        int(rot))
                        l_ol.append([l_olx, l_oly])

                str_ol = ''
                c = 0
                for i in range(0, len(n_im)):
                    if len(n_im[i]) == 1:
                        str_ol = str_ol + 'For S' + str(i + 1) .zfill(3) + ', {:.1f}'.format(l_ol[i]) + '%\n'
                        for j in range(0, n_im[i][0]):
                            pl.write('S' + str(i + 1) .zfill(3) + 'I' + str(j + 1).zfill(4)
                                     + ',{:.32f},{:.32f},{:.32f}'.format(v_x_o[c], v_y_o[c], v_z_o[c])
                                     + ',{:.32f},{:.32f}'.format(v_t[0], v_r[0])
                                     + ',0.00000000000000000000000000000000,'
                                     + '{:.6f},{:.8f}'.format(mag[i], v_wd_o[c])
                                     + '\n')
                            c = c + 1
                    elif len(n_im[i]) == 2:
                        str_ol = str_ol + 'For S' + str(i + 1).zfill(3) + ', {:.1f}'.format(l_ol[i][0]) + '% in x\n' + \
                                 'For S' + str(i + 1).zfill(3) + ', {:.1f}'.format(l_ol[i][1]) + '% in y\n'
                        for j in range(0, n_im[i][0]*n_im[i][1]):
                            pl.write('S' + str(i + 1).zfill(3) + 'I' + str(j + 1).zfill(4)
                                     + ',{:.32f},{:.32f},{:.32f}'.format(v_x_o[c], v_y_o[c], v_z_o[c])
                                     + ',{:.32f},{:.32f}'.format(v_t[0], v_r[0])
                                     + ',0.00000000000000000000000000000000,'
                                     + '{:.6f},{:.8f}'.format(mag[i], v_wd_o[c])
                                     + '\n')
                            c = c + 1

                    # if i > 0:
                    #     ax1.arrow(v_x_f[i-1], v_y_f[i-1], v_x_i[i]-v_x_f[i-1], v_y_i[i]-v_y_f[i-1],
                    #               head_width=2, head_length=3, fc='k', ec='k', length_includes_head=True)

                pl.write('S' + str(totim + 1).zfill(4) + 'end'
                         + ',{:.32f},{:.32f},{:.32f}'.format(v_x_o[-1], v_y_o[-1], v_z_o[-1])
                         + ',{:.32f},{:.32f}'.format(v_t[-1], v_r[-1])
                         + ',0.00000000000000000000000000000000,'
                         + '{:.6f},{:.8f}'.format(mag[-1], v_wd_o[-1])
                         + '\n')

            ax1.set_xlim([25, 105])
            ax1.set_ylim([25, 105])
            ax1.set_aspect('equal')
            ax1.set_xlabel('X [mm]')
            ax1.set_ylabel('Y [mm]')

            fig_c_a = draw_figure(window['figCanvas'].TKCanvas, fig)

            if event in 'Generate point list':
                window['time'].update('{:.2f}'.format(totim*0.5) +
                                      ' | ' + '{:.2f}'.format(totim*2) +
                                      ' minutes if "slow" | "very slow"')
                if values['ol_cb']:
                    ol_window = open_ol_window(len(n_im))
                    ol_window['output'].update(str_ol)

                    cnim = 0
                    for i in range(0, len(n_im)):
                        if len(n_im[i]) == 1:
                            cnim_end = cnim+n_im[i][0]
                        elif len(n_im[i]) == 2:
                            cnim_end = cnim+n_im[i][0]*n_im[i][1]
                        else:
                            cnim_end = 0

                        olf = plt.figure()
                        oax = olf.gca()
                        oax.plot(1000 * v_x_o[cnim:cnim_end], 1000 * v_y_o[cnim:cnim_end], 'x-')
                        draw_rectangle(v_x_o, v_y_o, cnim, cnim_end, int(rot), mag[i], oax)

                        oax.set_aspect('equal')
                        oax.set_xlabel('X [mm]')
                        oax.set_ylabel('Y [mm]')
                        draw_figure(ol_window['olCanvas'+str(i)].TKCanvas, olf)

                        cnim = cnim_end

        else:
            sg.Popup('Input error', 'The # of images and the magnification of each series must be entered separated by '
                                    'a comma')
    else:
        sg.Popup('Input error', 'Tilt or rotation of first and last point of a serie are not equal.\n'
                                'Please enter valid initial points')
window.close()
