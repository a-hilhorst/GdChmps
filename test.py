import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import PySimpleGUI as Sg
import os

ptListIn = os.getcwd() + '\IN_2_4_2.TXT'
ptListOut = os.getcwd() + '\IN_2_4_2_OUT.TXT'

layout = [
    [Sg.Text('# images of each series: ',
             tooltip='HELLO',
             size=(20, 1)), Sg.Input(key='n_im')],
    [Sg.Text('Magnification of each series: ', tooltip='HELLO',
             size=(20, 1)), Sg.Input(key='mag')],
    [Sg.Submit('T-E-S-T'), Sg.Exit()]
]

window = Sg.Window('Test', layout)

while True:
    event, values = window.read()
    if event in (None, 'Exit'):
        break

    n_im = []
    mag = []
    while True:
        try:
            n_im = values['n_im'].split(',')
            n_im = [i.split('x') for i in n_im]
            n_im = [[int(i) for i in j] for j in n_im]
            mag = [float(i) for i in values['mag'].split(',')]
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
        print('OK')
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

                    for k in range(lv_x_o1, lv_x_o2):
                        rect = patches.Rectangle((1000 * v_x_o[k], 1000 * v_y_o[k]), 114/mag[i], 85.7/mag[i], linewidth=1,
                                                 edgecolor='k', facecolor='none', alpha=0.5)
                        ax1.add_patch(rect)

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

                    lv_x_o2 = len(v_x_o)
                    ax1.plot(1000 * v_x_o[lv_x_o1:lv_x_o2], 1000 * v_y_o[lv_x_o1:lv_x_o2], 'x-')
                    c = c + 4
                else:
                    Sg.Popup('not 1x1')

        ax1.set_xlim([25, 105])
        ax1.set_ylim([25, 105])
        ax1.set_aspect('equal')
        ax1.set_xlabel('X [mm]')
        ax1.set_ylabel('Y [mm]')

        plt.show(block=False)

    else:
        Sg.Popup(str(len(v_x)) + ' ' + str(2*lnim))

    if event in 'T-E-S-T':
        Sg.Popup(n_im)

