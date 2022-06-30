import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as transforms
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as Sg
import os

ptListIn = os.getcwd() + '/IN_2_4_2.TXT'
ptListOut = os.getcwd() + '/IN_2_4_2_OUT.TXT'


def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg


# vC_x = 114 / mag[i] / 2
# vC_y = 85.7 / mag[i] / 2
# vC_x_p = vC_x * np.cos(int(rot) * np.pi / 180) - vC_y * np.sin(int(rot) * np.pi / 180)
# vC_y_p = vC_x * np.sin(int(rot) * np.pi / 180) + vC_y * np.cos(int(rot) * np.pi / 180)
# t_t = transforms.Affine2D().translate(1000 * v_x_o[lv_x_o1 + 2*n_im[i][0]-1] - vC_x_p,
#                                       1000 * v_y_o[lv_x_o1 + 2*n_im[i][0]-1] - vC_y_p)


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

    ol = 100*si/su  # overlap in %
    return ol


column1 = [
    [Sg.Text('# images of each series: ',
             tooltip='HELLO',
             size=(20, 1)), Sg.Input(key='n_im', default_text='10,2x5,5')],
    [Sg.Text('Magnification of each series: ', tooltip='HELLO',
             size=(20, 1)), Sg.Input(key='mag', default_text='10,20,20')],
    [Sg.Text('Scan rot :', tooltip='HELLO',
             size=(20, 1)), Sg.Input(key='rot', default_text='0')],
    [Sg.Submit('Submit')],
    [Sg.Text('Overlap :\n==============')],
    [Sg.Text('',
             size=(20, 1), key='output')],
    [Sg.Exit()]
]
column2 = [
    [Sg.Canvas(key='figCanvas')]
]
layout = [
    [Sg.Column(column1),
     Sg.VSeparator(),
     Sg.Column(column2)]
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

window = Sg.Window('Test', layout, finalize=True, resizable=True, element_justification="left")
fig_c_a = draw_figure(window['figCanvas'].TKCanvas, fig)

while True:
    event, values = window.read()
    if event in (None, 'Exit'):
        break

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
            Sg.Popup('Error', 'Please enter # of images as integers separated by a comma\n'
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
            for i in range(0, len(n_im)):
                if len(n_im[i]) == 1:
                    # line A to B
                    lv_x_o1 = len(v_x_o)

                    v_x_o = np.concatenate((v_x_o, np.linspace(v_x[c], v_x[c+1], n_im[i][0])), axis=0)
                    v_y_o = np.concatenate((v_y_o, np.linspace(v_y[c], v_y[c+1], n_im[i][0])), axis=0)
                    v_z_o = np.concatenate((v_z_o, np.linspace(v_z[c], v_z[c+1], n_im[i][0])), axis=0)
                    v_wd_o = np.concatenate((v_wd_o, np.linspace(v_wd[c], v_wd[c+1], n_im[i][0])), axis=0)
                    c = c + 2

                    lv_x_o2 = len(v_x_o)
                    ax1.plot(1000 * v_x_o[lv_x_o1:lv_x_o2], 1000 * v_y_o[lv_x_o1:lv_x_o2], 'x-')

                    t_r = transforms.Affine2D().rotate_deg(int(rot))
                    for k in range(lv_x_o1, lv_x_o2):
                        rect = patches.Rectangle((0, 0), 114 / mag[i], 85.7 / mag[i],
                                                 linewidth=1,
                                                 edgecolor='k',
                                                 facecolor='none',
                                                 alpha=0.5)

                        vC_x = 114 / mag[i] / 2
                        vC_y = 85.7 / mag[i] / 2
                        vC_x_p = vC_x * np.cos(int(rot) * np.pi / 180) - vC_y * np.sin(int(rot) * np.pi / 180)
                        vC_y_p = vC_x * np.sin(int(rot) * np.pi / 180) + vC_y * np.cos(int(rot) * np.pi / 180)
                        t_t = transforms.Affine2D().translate(1000 * v_x_o[k] - vC_x_p,
                                                              1000 * v_y_o[k] - vC_y_p)
                        rect.set_transform(t_r + t_t + ax1.transData)

                        ax1.add_patch(rect)

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

                    t_r = transforms.Affine2D().rotate_deg(int(rot))
                    for k in range(lv_x_o1, lv_x_o2):
                        rect = patches.Rectangle((0, 0), 114 / mag[i], 85.7 / mag[i],
                                                 linewidth=1,
                                                 edgecolor='k',
                                                 facecolor='none',
                                                 alpha=0.5)

                        vC_x = 114 / mag[i] / 2
                        vC_y = 85.7 / mag[i] / 2
                        vC_x_p = vC_x * np.cos(int(rot) * np.pi / 180) - vC_y * np.sin(int(rot) * np.pi / 180)
                        vC_y_p = vC_x * np.sin(int(rot) * np.pi / 180) + vC_y * np.cos(int(rot) * np.pi / 180)
                        t_t = transforms.Affine2D().translate(1000 * v_x_o[k] - vC_x_p,
                                                              1000 * v_y_o[k] - vC_y_p)
                        rect.set_transform(t_r + t_t + ax1.transData)

                        ax1.add_patch(rect)

                    l_olx = overlap(mag[i],
                                    [1000 * (v_x_o[lv_x_o1] - v_x_o[lv_x_o1 + 2 * n_im[i][0] - 1]),
                                     1000 * (v_y_o[lv_x_o1] - v_y_o[lv_x_o1 + 2 * n_im[i][0] - 1])],
                                    int(rot))
                    l_oly = overlap(mag[i],
                                    [1000 * (v_x_o[lv_x_o1] - v_x_o[lv_x_o1 + 1]),
                                     1000 * (v_y_o[lv_x_o1] - v_y_o[lv_x_o1 + 1])],
                                    int(rot))
                    l_ol.append([l_olx, l_oly])

                    #  r
                    rect = patches.Rectangle((0, 0), 114 / mag[i], 85.7 / mag[i],
                                             linewidth=2,
                                             edgecolor='r',
                                             facecolor='none',
                                             alpha=1)

                    vC_x = 114 / mag[i] / 2
                    vC_y = 85.7 / mag[i] / 2
                    vC_x_p = vC_x * np.cos(int(rot) * np.pi / 180) - vC_y * np.sin(int(rot) * np.pi / 180)
                    vC_y_p = vC_x * np.sin(int(rot) * np.pi / 180) + vC_y * np.cos(int(rot) * np.pi / 180)
                    t_t = transforms.Affine2D().translate(1000 * v_x_o[lv_x_o1 + 2*n_im[i][0]-1] - vC_x_p,
                                                          1000 * v_y_o[lv_x_o1 + 2*n_im[i][0]-1] - vC_y_p)
                    rect.set_transform(t_r + t_t + ax1.transData)

                    ax1.add_patch(rect)

                    #  b
                    rect = patches.Rectangle((0, 0), 114 / mag[i], 85.7 / mag[i],
                                             linewidth=2,
                                             edgecolor='b',
                                             facecolor='none',
                                             alpha=1)

                    vC_x = 114 / mag[i] / 2
                    vC_y = 85.7 / mag[i] / 2
                    vC_x_p = vC_x * np.cos(int(rot) * np.pi / 180) - vC_y * np.sin(int(rot) * np.pi / 180)
                    vC_y_p = vC_x * np.sin(int(rot) * np.pi / 180) + vC_y * np.cos(int(rot) * np.pi / 180)
                    t_t = transforms.Affine2D().translate(1000 * v_x_o[lv_x_o1] - vC_x_p,
                                                          1000 * v_y_o[lv_x_o1] - vC_y_p)
                    rect.set_transform(t_r + t_t + ax1.transData)

                    ax1.add_patch(rect)

                    #  g
                    rect = patches.Rectangle((0, 0), 114 / mag[i], 85.7 / mag[i],
                                             linewidth=2,
                                             edgecolor='g',
                                             facecolor='none',
                                             alpha=1)

                    vC_x = 114 / mag[i] / 2
                    vC_y = 85.7 / mag[i] / 2
                    vC_x_p = vC_x * np.cos(int(rot) * np.pi / 180) - vC_y * np.sin(int(rot) * np.pi / 180)
                    vC_y_p = vC_x * np.sin(int(rot) * np.pi / 180) + vC_y * np.cos(int(rot) * np.pi / 180)
                    t_t = transforms.Affine2D().translate(1000 * v_x_o[lv_x_o1 + 1] - vC_x_p,
                                                          1000 * v_y_o[lv_x_o1 + 1] - vC_y_p)
                    rect.set_transform(t_r + t_t + ax1.transData)

                    ax1.add_patch(rect)

                else:
                    Sg.Popup('not 1x1')

        # plt.show(block=False)
        ax1.set_xlim([25, 105])
        ax1.set_ylim([25, 105])
        ax1.set_aspect('equal')
        ax1.set_xlabel('X [mm]')
        ax1.set_ylabel('Y [mm]')

        fig_c_a = draw_figure(window['figCanvas'].TKCanvas, fig)

    else:
        Sg.Popup(str(len(v_x)) + ' ' + str(2*lnim))

    if event in 'Submit':
        # Sg.Popup(n_im)
        window['output'].update(str(len(v_x_o)))
        Sg.Popup(l_ol)

