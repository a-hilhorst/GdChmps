import numpy as np
# import matplotlib.pyplot as plt
import PySimpleGUI as Sg
import os

ptListIn = os.getcwd() + '\IN_2_4_2.TXT'
ptListOut = os.getcwd() + '\IN_2_4_2_OUT.TXT'

layout = [
    [Sg.Text('# images of each series: ',
             tooltip='Enter the number of images of each series as an integer separated by a comma \n'
                     'e.g. 2,3,4 for 3 series of 2, 3, and 4 images respectively.',
             size=(20, 1)), Sg.Input(key='n_im')],
    [Sg.Submit('T-E-S-T'), Sg.Exit()]
]

window = Sg.Window('Test', layout)

while True:
    event, values = window.read()
    if event in (None, 'Exit'):
        break

    n_im = []
    while True:
        try:
            n_im = values['n_im'].split(',')
            n_im = [i.split('x') for i in n_im]
            n_im = [[int(i) for i in j] for j in n_im]
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

        with open(ptListOut, 'w') as pl:
            pl.write('Leo Points List\n')
            pl.write('Absolute\n')
            pl.write('Label,X,Y,Z,T,R,M,Mag,WD\n')
            pl.write(str(totim + 1) + '\n')
    else:
        Sg.Popup(str(len(v_x)) + ' ' + str(2*lnim))

    if event in 'T-E-S-T':
        Sg.Popup(n_im)

