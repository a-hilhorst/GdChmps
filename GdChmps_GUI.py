import numpy as np
import matplotlib.pyplot as plt
import PySimpleGUI as Sg

layout = [
    [Sg.Text('Scan Rotate is not taken into account.')],
    [Sg.Text('Set beam shift to 0 before taking initial and end points.')],
    [Sg.Text('Tilt and Rotation must not change between the initial and final points of a serie.')],
    [Sg.Text('Initial points:', tooltip='Point Lists are saved in C:\\ProgramData\\Carl Zeiss\\SmartSEM\\User\\sem',
             size=(20, 1))],
    [Sg.Input(key='ptInit', size=(60, 1)), Sg.FileBrowse(target='ptInit', size=(6, 1),
                                                         initial_folder='C:/ProgramData/Carl Zeiss/SmartSEM/User/sem')],
    [Sg.Text('# images of each series: ',
             tooltip='Enter the number of images of each series as an integer separated by a comma\n'
                     'e.g. 2,3,4 for 3 series of 2, 3, and 4 images respectively.',
             size=(20, 1)), Sg.Input(key='n_im')],
    [Sg.Text('Magnification of each series: ', tooltip='Enter the magnification of each series separated by a comma\n'
                                                       'e.g. 1000,3000,2000 for 3 series with mag = 1k, 3k, and 2k '
                                                       'respectively.',
             size=(20, 1)), Sg.Input(key='mag')],
    [Sg.Text('Save as: ', size=(20, 1))],
    [Sg.Input(key='saveName', size=(58, 1)), Sg.SaveAs(target='saveName', size=(8, 1))],
    [Sg.Text('In Stage Points List > On Goto, check Move XY Only.')],
    [Sg.Submit('Generate point list'), Sg.Exit()]
]

window = Sg.Window('Point List Generator', layout)

while True:
    event, values = window.read()
    if event in (None, 'Exit'):
        break

    ptListIn = values['ptInit']
    ptListOut = values['saveName']

    n_im = []
    mag = []
    while True:
        try:
            n_im = [int(i) for i in values['n_im'].split(',')]
            mag = [float(i) for i in values['mag'].split(',')]
        except ValueError:
            Sg.Popup('Error', 'Please enter # of images as integers separated by a comma\n'
                              'Please enter magnification as doubles separated by a comma')
            event, values = window.read()
            continue
        break

    with open(ptListIn, 'r') as pl:
        pll = pl.readlines()

        v_x_i = []
        v_y_i = []
        v_z_i = []
        v_t_i = []
        v_r_i = []
        v_wd_i = []
        v_x_f = []
        v_y_f = []
        v_z_f = []
        v_t_f = []
        v_r_f = []
        v_wd_f = []
        for i in range(4, len(pll)):
            tmp = pll[i].split(',')
            if (i % 2) == 0:
                v_x_i.append(float(tmp[1]))
                v_y_i.append(float(tmp[2]))
                v_z_i.append(float(tmp[3]))
                v_t_i.append(float(tmp[4]))
                v_r_i.append(float(tmp[5]))
                v_wd_i.append(float(tmp[8][0:-1]))
            else:
                v_x_f.append(float(tmp[1]))
                v_y_f.append(float(tmp[2]))
                v_z_f.append(float(tmp[3]))
                v_t_f.append(float(tmp[4]))
                v_r_f.append(float(tmp[5]))
                v_wd_f.append(float(tmp[8][0:-1]))

    deltax = ''
    deltay = ''
    deltaz = ''
    if (v_t_i == v_t_f) & (v_r_i == v_r_f):
        if (len(n_im) == len(mag)) & (len(n_im) == len(v_x_f)):
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
                pl.write(str(sum(n_im) + 1) + '\n')

                for i in range(0, len(v_x_i)):
                    v_x = np.linspace(v_x_i[i], v_x_f[i], n_im[i])
                    v_y = np.linspace(v_y_i[i], v_y_f[i], n_im[i])
                    v_z = np.linspace(v_z_i[i], v_z_f[i], n_im[i])
                    v_wd = np.linspace(v_wd_i[i], v_wd_f[i], n_im[i])

                    ax1.plot(1000*v_x, 1000*v_y, 'x-')
                    # ax2.plot(1000*v_x, 1000*v_z, 'x-')

                    if len(v_x) > 1:
                        deltax += ', {:.2f} \u03BC'.format(1000000*(v_x[1] - v_x[0])) + 'm'
                        deltay += ', {:.2f} \u03BC'.format(1000000*(v_y[1] - v_y[0])) + 'm'
                        deltaz += ', {:.2f} \u03BC'.format(1000000*(v_z[1] - v_z[0])) + 'm'

                    for j in range(0, n_im[i]):
                        if n_im[i] > 99:
                            if j > 98:
                                pl.write('S' + str(i + 1) + 'I' + str(j + 1)
                                         + ',{:.32f},{:.32f},{:.32f}'.format(v_x[j], v_y[j], v_z[j])
                                         + ',{:.32f},{:.32f}'.format(v_t_i[i], v_r_i[i])
                                         + ',0.00000000000000000000000000000000,'
                                         + '{:.6f},{:.8f}'.format(mag[i], v_wd[j])
                                         + '\n')
                            elif j > 8:
                                pl.write('S' + str(i + 1) + 'I0' + str(j + 1)
                                         + ',{:.32f},{:.32f},{:.32f}'.format(v_x[j], v_y[j], v_z[j])
                                         + ',{:.32f},{:.32f}'.format(v_t_i[i], v_r_i[i])
                                         + ',0.00000000000000000000000000000000,'
                                         + '{:.6f},{:.8f}'.format(mag[i], v_wd[j])
                                         + '\n')
                            else:
                                pl.write('S' + str(i + 1) + 'I00' + str(j + 1)
                                         + ',{:.32f},{:.32f},{:.32f}'.format(v_x[j], v_y[j], v_z[j])
                                         + ',{:.32f},{:.32f}'.format(v_t_i[i], v_r_i[i])
                                         + ',0.00000000000000000000000000000000,'
                                         + '{:.6f},{:.8f}'.format(mag[i], v_wd[j])
                                         + '\n')
                        elif n_im[i] > 9:
                            if j > 8:
                                pl.write('S' + str(i + 1) + 'I' + str(j + 1)
                                         + ',{:.32f},{:.32f},{:.32f}'.format(v_x[j], v_y[j], v_z[j])
                                         + ',{:.32f},{:.32f}'.format(v_t_i[i], v_r_i[i])
                                         + ',0.00000000000000000000000000000000,'
                                         + '{:.6f},{:.8f}'.format(mag[i], v_wd[j])
                                         + '\n')
                            else:
                                pl.write('S' + str(i + 1) + 'I0' + str(j + 1)
                                         + ',{:.32f},{:.32f},{:.32f}'.format(v_x[j], v_y[j], v_z[j])
                                         + ',{:.32f},{:.32f}'.format(v_t_i[i], v_r_i[i])
                                         + ',0.00000000000000000000000000000000,'
                                         + '{:.6f},{:.8f}'.format(mag[i], v_wd[j])
                                         + '\n')
                        else:
                            pl.write('S' + str(i + 1) + 'I' + str(j + 1)
                                     + ',{:.32f},{:.32f},{:.32f}'.format(v_x[j], v_y[j], v_z[j])
                                     + ',{:.32f},{:.32f}'.format(v_t_i[i], v_r_i[i])
                                     + ',0.00000000000000000000000000000000,'
                                     + '{:.6f},{:.8f}'.format(mag[i], v_wd[j])
                                     + '\n')

                    if i > 0:
                        ax1.arrow(v_x_f[i-1], v_y_f[i-1], v_x_i[i]-v_x_f[i-1], v_y_i[i]-v_y_f[i-1],
                                  head_width=2, head_length=3, fc='k', ec='k', length_includes_head=True)

                pl.write('S' + str(len(v_x_i) + 1) + 'end'
                         + ',{:.32f},{:.32f},{:.32f}'.format(v_x_f[-1], v_y_f[-1], v_z_f[-1])
                         + ',{:.32f},{:.32f}'.format(v_t_i[-1], v_r_i[-1])
                         + ',0.00000000000000000000000000000000,'
                         + '{:.6f},{:.8f}'.format(mag[-1], v_wd_f[-1])
                         + '\n')

            ax1.set_xlim([25, 105])
            ax1.set_ylim([25, 105])
            ax1.set_aspect('equal')
            ax1.set_xlabel('X [mm]')
            ax1.set_ylabel('Y [mm]')
            # ax1.invert_xaxis()

            plt.show(block=False)

            if len(v_x) > 1:
                if event in 'Generate point list':
                    Sg.Popup('Displacements',
                             'Delta X for each series is:' + deltax[1:] + '\n' +
                             'Delta Y for each series is:' + deltay[1:] + '\n' +
                             'Delta Z for each series is:' + deltaz[1:])
        else:
            Sg.Popup('Input error', 'The # of images and the magnification of each series must be entered separated by '
                                    'a comma')
    else:
        Sg.Popup('Input error', 'Tilt or rotation of first and last point of a serie are not equal.\n'
                                'Please enter valid initial points')
window.close()
