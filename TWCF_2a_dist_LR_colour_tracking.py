#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Distance comparison across blind spot
TWCF IIT vs PP experiment 2a piloting
Authors: Clement Abbatecola, Belén María Montabes de la Cruz
    Code version:
        0.7 # 2022/10/25 fixation before start of the first run ('up' key to start)
        0.6 # 2022/09/22 introducing system to read bs and color parameters from files
        0.5 # 2022/09/18 bug fixes + added abort trial option
        0.4 # 2022/09/06 introducing colour/monocular stimulation
        0.3 # 2022/08/30 simplifying, debugging
        0.2 # final version for piloting
"""

from psychopy import prefs
prefs.hardware['audioLib'] = ['PTB'] # not used actually

from psychopy import core, visual, gui, data, event, monitors, sound
from psychopy.tools.coordinatetools import pol2cart, cart2pol
import numpy as np
import random, datetime, os
from glob import glob
from itertools import compress

# custom sources:
from fusion_stim import fusionStim
from LT import *

# for method of constant stimuli:
# (but we decided on staircases)

# def getStimuli(hemisphere='left'):

#     eyes = ['left','righ']
#     pos_idc = [0,1,2,3]
#     intervals = [-3.5, -1.5, -0.5, 0.5, 1.5, 3.5]

#     # interval_repetitions = [5, 5, 10, 10, 5, 5]

#     stimuli = []

#     for eye in eyes:
#         if eye == hemisphere:
#             interval_repetitions = [4, 4, 8, 8, 4, 4]
#         else:
#             interval_repetitions = [4, 4, 8, 8, 4, 4]
#         for pos_idx in pos_idc:
#             for ino in range(len(intervals)):
#                 interval = intervals[ino]
#                 repetitions = interval_repetitions[ino]

#                 stim = {'eye':eye, 'pos':pos_idx, 'interval':interval}
#                 stimuli += [stim]*repetitions

#     return(stimuli)

######
#### Initialize experiment
######

## parameters
nRevs   = 10   #
nTrials  = 30  # might need to tweak these values????
letter_height = 1

## path
main_path = "."
data_path = main_path + "/data/"

# get experiment info?
expInfo = {'ID':'test',  'hemifield':['left','right'], 'track eyes':['both','left','right','off']}
dlg = gui.DlgFromDict(expInfo, title='Infos')

hemifield = expInfo['hemifield']

track_eyes = expInfo['track eyes']

trackLeftEye, trackRightEye = {'both':  [True,  True ],
                               'left':  [True,  False],
                               'right': [False, True ],
                               'off':   [False, False] }[track_eyes]

random.seed(expInfo['ID'] + hemifield) # random sequences unique to each participant and run


## files
run_no = 1 #write non-existent output file
filename = 'dist_' + hemifield + '_' + expInfo['ID'].lower() + '_'
while (filename + str(run_no) + '.txt') in os.listdir(data_path): run_no += 1
respFile = open(data_path + filename + str(run_no) + '.txt','w')

if any([trackLeftEye, trackRightEye]):
    # make eye-tracking directory for participant:
    eyetrack_path = '%seyetracking/%s_%s_%d/'%(data_path, expInfo['ID'], hemifield, run_no)
    if not os.path.exists(eyetrack_path):
        os.makedirs(eyetrack_path)

## blindspot parameters
bs_file = open(glob(main_path + "/mapping_data/" + expInfo['ID'].lower() + {'left':"_LH_blindspot*.txt", 'right': "_RH_blindspot*.txt"}[hemifield])[-1],'r')
bs_param = bs_file.read().replace('\t','\n').split('\n')
bs_file.close()
if hemifield == 'left':
    spot_left_cart = eval(bs_param[1])
    spot_left = cart2pol(spot_left_cart[0], spot_left_cart[1])
    spot_left_size = eval(bs_param[3])
if hemifield == 'right':
    spot_right_cart = eval(bs_param[1])
    spot_right = cart2pol(spot_right_cart[0], spot_right_cart[1])
    spot_right_size = eval(bs_param[3])

# scale blindspot marker
# we want half the surface?

bs_scale = 0.5
bsf = bs_scale**0.5
if hemifield == 'left':
    marker_size = [spot_left_size[0] * bsf, spot_left_size[1] * bsf]
if hemifield == 'right':
    marker_size = [spot_right_size[0] * bsf, spot_right_size[1] * bsf]


'''
distance of reference between dots (target)
=> width of blindspot + 2 (dot width, padding) + 2 (to account for a max jitter of 1 on either side)
'''
# plus 1 more dva for padding (marius)
if hemifield == 'left':
    tar =  spot_left_size[0] + 2 + 2 + 1
    # angle between center and top/bottom of blindspot + 4 (dot width, padding)
    if spot_left_cart[1] < 0:
        ang_up = (cart2pol(spot_left_cart[0], spot_left_cart[1] - spot_left_size[1])[0] - spot_left[0])
    else:
       ang_up = (cart2pol(spot_left_cart[0], spot_left_cart[1] + spot_left_size[1])[0] - spot_left[0])
if hemifield == 'right':
    tar =  spot_right_size[0] + 2 + 2 + 1
    # angle between center and top/bottom of blindspot + 4 (dot width, padding)
    if spot_right_cart[1] < 0:
        ang_up = (cart2pol(spot_right_cart[0], spot_right_cart[1] - spot_right_size[1])[0] - spot_right[0])
    else:
       ang_up = (cart2pol(spot_right_cart[0], spot_right_cart[1] + spot_right_size[1])[0] - spot_right[0])



## colour (eye) parameters
col_file = open(glob(main_path + "/mapping_data/" + expInfo['ID'] + "_col_cal*.txt")[-1],'r')
col_param = col_file.read().replace('\t','\n').split('\n')
col_file.close()
col_back  = eval(col_param[1])
col_left  = eval(col_param[3])
col_right = eval(col_param[5])
col_both  = [col_left[0], -1, col_right[2]]



gammaGrid = np.array([  [  0., 135.44739,  2.4203537, np.nan, np.nan, np.nan  ],
                        [  0.,  27.722954, 2.4203537, np.nan, np.nan, np.nan  ],
                        [  0.,  97.999275, 2.4203537, np.nan, np.nan, np.nan  ],
                        [  0.,   9.235623, 2.4203537, np.nan, np.nan, np.nan  ]], dtype=float)

resolution = [1920, 1080]
size = [59.8, 33.6]
distance = 50


mymonitor = monitors.Monitor(name='temp',
                             distance=distance,
                             width=size[0])
mymonitor.setGammaGrid(gammaGrid)
mymonitor.setSizePix(resolution)




## window & elements
# win = visual.Window([2500,1000],allowGUI=True, monitor='ccni', units='deg', viewPos = [0,0], fullscr = True, color=[.5,-1,.5])
# win = visual.Window([1500,800],allowGUI=True, monitor='testMonitor', units='deg', viewPos = [13,0], fullscr = False, color=[.5,-1,.5])

# win = visual.Window(resolution,allowGUI=True, monitor=mymonitor, units='deg', viewPos = [13,0], fullscr = False, color=[.5,-1,.5])
win = visual.Window(resolution, allowGUI=True, monitor=mymonitor, units='deg', fullscr = True, color=col_back, screen=1) # back to same lay-out as in blindspot mapping task, to keep stuff really aligned

win.mouseVisible = False
fixation = visual.ShapeStim(win, 
                            vertices = ((-1., -1.), (1., 1.), (0,0), (-1., 1.), (1., -1.)), 
                            lineWidth = 2, 
                            units = 'deg', 
                            size = (0.36, 0.36), 
                            closeShape = False, 
                            lineColor = 'black')

class myHash(win,
             pos = [0,0],
             units='norm',
             size=(0.1,0.1),
             lineWidth = 2,
             lineColor='black',
             ori=15):

    def __init__(self):
        self.win = win
        self.pos = pos
        self.units = units
        self.size = size
        self.lineWidth = lineWidth
        self.lineColor = lineColor
        self.ori = ori

        self.setLines()

    def setLines(self)
        self.lines = []
        ex = self.size[0]/2
        ey = self.size[1]/2
        self.lines.append( visual.Line(win   = self.win,
                                       start = [ex*(-1), ey*(-1/3)],
                                       end   = [ex*( 1), ey*(-1/3)],
                                       units = self.units, pos = self.pos, lineWidth = self.lineWidth, lineColor = self.lineColor, ori = self.ori) )
        self.lines.append( visual.Line(win   = self.win,
                                       start = [ex*(-1), ey*( 1/3)],
                                       end   = [ex*( 1), ey*( 1/3)],
                                       units = self.units, pos = self.pos, lineWidth = self.lineWidth, lineColor = self.lineColor, ori = self.ori) )
        self.lines.append( visual.Line(win   = self.win,
                                       start = [ex*(-1/3), ey*(-1)],
                                       end   = [ex*(-1/3), ey*( 1)],
                                       units = self.units, pos = self.pos, lineWidth = self.lineWidth, lineColor = self.lineColor, ori = self.ori) )
        self.lines.append( visual.Line(win   = self.win,
                                       start = [ex*( 1/3), ey*(-1)],
                                       end   = [ex*( 1/3), ey*( 1)],
                                       units = self.units, pos = self.pos, lineWidth = self.lineWidth, lineColor = self.lineColor, ori = self.ori) )

    def draw(self):
        for (ln in self.lines):
            ln.draw()

hash_fix = myHash(win=win)

## instructions
visual.TextStim(win,'Troughout the experiment you will fixate at a cross that will be located at the middle of the screen.   \
It is important that you fixate on this cross at all times.\n\n You will be presented with pairs of dots. You will have to indicate which dots were closer together.\n\n Left arrow = first pair of dots were closer together.\
\n\n Right arrow = second pair of dots were closer together.\n\n\n Press the space bar to start the experiment.', height = letter_height,wrapWidth=30, color = 'black').draw()
win.flip()
k = ['wait']
while k[0] not in ['escape','space']:
    k = event.waitKeys()
if k[0] in ['escape']:
    win.close()
    core.quit()


cfg = {}
cfg['hw'] = {}
cfg['hw']['win'] = win

######
#### Prepare stimulation
######


# tick = sound.Sound('268108__nenadsimic__button-tick.wav', sampleRate=44100)


## stimuli
point_1 = visual.Circle(win, radius = .5, pos = pol2cart( 0, 3), fillColor = col_both, lineColor = None)
point_2 = visual.Circle(win, radius = .5, pos = pol2cart( 0, 6), fillColor = col_both, lineColor = None)
point_3 = visual.Circle(win, radius = .5, pos = pol2cart(45, 3), fillColor = col_both, lineColor = None)
point_4 = visual.Circle(win, radius = .5, pos = pol2cart(45, 6), fillColor = col_both, lineColor = None)


# fusion patches:
hiFusion = fusionStim(win=win, units='deg', pos=[0, 6], square=0.3)
loFusion = fusionStim(win=win, units='deg', pos=[0,-6], square=0.3)

if hemifield == 'left':
    blindspot = visual.Circle(win, radius = .5, pos = [7,0], fillColor=col_left, lineColor = None)
    blindspot.pos = spot_left_cart
if hemifield == 'right':
    blindspot = visual.Circle(win, radius = .5, pos = [7,0], fillColor=col_right, lineColor = None)
    blindspot.pos = spot_right_cart
blindspot.size = marker_size 
# blindspot.size = spot_left_size

## prepare trials
if hemifield == 'left':
    positions = {
        "left-top": [(spot_left[0] - ang_up, spot_left[1] - tar/2), (spot_left[0] - ang_up, spot_left[1] + tar/2)],
        "left-mid": [(spot_left[0] +     00, spot_left[1] - tar/2), (spot_left[0] +     00, spot_left[1] + tar/2)],
        "left-bot": [(spot_left[0] + ang_up, spot_left[1] - tar/2), (spot_left[0] + ang_up, spot_left[1] + tar/2)],
    }
    # First column is target, second column is foil, those are the valid trials in the proportion we want
    pos_array = [["left-mid", "left-top"],
                 ["left-mid", "left-bot"],
                 ["left-top", "left-bot"],
                 ["left-bot", "left-top"]]

if hemifield == 'right':
    positions = {
        "right-top": [(spot_right[0] + ang_up, spot_right[1] - tar/2), (spot_right[0] + ang_up, spot_right[1] + tar/2)],
        "right-mid": [(spot_right[0] +     00, spot_right[1] - tar/2), (spot_right[0] +     00, spot_right[1] + tar/2)],
        "right-bot": [(spot_right[0] - ang_up, spot_right[1] - tar/2), (spot_right[0] - ang_up, spot_right[1] + tar/2)],
    }
    # First column is target, second column is foil, those are the valid trials in the proportion we want
    pos_array = [["right-mid", "right-top"],
                 ["right-mid", "right-bot"],
                 ["right-top", "right-bot"],
                 ["right-bot", "right-top"]]

random.shuffle(pos_array)

respFile.write(''.join(map(str, ["Start: \t" + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "\n"])))
respFile.write("Block\tTrial\tResp\tTarg_loc\tFoil_loc\tTarg_len\tDifference\tWhich_first\tCorrect\tReversal\tFoil_type\tEye\tStair\n")
print(datetime.datetime.now().strftime("%Y-%m-%d-%H-%M"))
print("Block", "Trial", "Resp", "Targ_loc", "Foil_loc", "Targ_len", "Difference", "Which_first", "Correct", "Reversal", "Foil_type", "Eye", "Stair")

# fixation help
fixation.draw()
win.flip()
k = ['wait']
while k[0] not in ['up', 'escape']:
    k = event.waitKeys()
if k[0] in ['escape']:
    win.close()
    core.quit()

# # # # # #
# Staircase
# # # # # #

intervals = [3.5, 3, 2.5, 2, 1.5, 1, .5, .1, -.1, -.5, -1, -1.5, -2, -2.5, -3, -3.5]
# this comes down to (target-foil) distance
# i.e. a high value means the target pair is longer than the foil pair

max_int_idx = len(intervals) - 1


meta_clock = core.Clock()
trial_clock = core.Clock()

trial = [0,0,0,0,0,0,0,0]
position = [pos_array[:],pos_array[:],pos_array[:],pos_array[:],pos_array[:],pos_array[:],pos_array[:],pos_array[:]]
foil_type = [1,-1,1,-1,1,-1,1,-1] # gets multiplied with interval difference... 
eye = ['left', 'left', 'right', 'right', 'left', 'left', 'right', 'right']
revs = [0,0,0,0,0,0,0,0]
direction = [1,1,1,1,-1,-1,-1,-1]
cur_int = [0,0,0,0,max_int_idx,max_int_idx,max_int_idx,max_int_idx] # current interval? if that is the index into the interval vector, they all start at 3.5 and none at -3.5
reversal = [False,False,False,False,False,False,False,False] # this is reset to a single boolean value later on, remove?
resps = [[True],[True],[True],[True],[True],[True],[True],[True]]
reversalIntensities = [[],[],[],[],[],[],[],[]] # this seems to be unused, remove?
stairs_ongoing = [True, True, True, True, True, True, True, True]

# for new staircases:
pairs_chosen = [[],[],[],[],[],[],[],[]]
last_step_directions = [0,0,0,0,0,0,0,0]




# marius_counter = 0

# stimuli = getStimuli()

# random.shuffle(stimuli)

abort = False




##### do eye-tracker calibration
if any([trackLeftEye, trackRightEye]):
    visual.TextStim(win,'SPACE starts eye-tracker calibration.', height = letter_height, color = 'white').draw()
    win.flip()

    waiting_for_response = True

    while waiting_for_response:
        keys = event.getKeys(keyList=['space'])
        if len(keys):
            if 'space' in keys:
                waiting_for_response = False
                win.flip()


    # initialise LiveTrack
    LiveTrack.Init()

    # Start LiveTrack raw data streaming (???)
    LiveTrack.SetResultsTypeRaw()

    # Start buffering data to the library
    LiveTrack.StartTracking()

    # left eye tracking only ???
    LiveTrack.SetTracking(leftEye=trackLeftEye,rightEye=trackRightEye)

    # do calibration:
    LTcal(cfg=cfg, trackLeftEye=trackLeftEye, trackRightEye=trackRightEye)

    # # set output to be calibrated output:
    # LiveTrack.SetResultsTypeCalibrated()

##### now start doing the task

visual.TextStim(win,'SPACE starts dot-distance task.', height = letter_height, color = 'white').draw()
win.flip()

waiting_for_response = True

while waiting_for_response:
    keys = event.getKeys(keyList=['space','escape'])
    if len(keys):
        if 'space' in keys:
            waiting_for_response = False

if any([trackLeftEye, trackRightEye]):
    # set output to be calibrated output:
    LiveTrack.SetResultsTypeCalibrated()

blocktrialcounter = 0
blockcounter = 1
trialcounter = 1

if any([trackLeftEye, trackRightEye]):

    calibrationcounter = 1

    LTstoreCalibration('%scalibration_%d.json'%(eyetrack_path, calibrationcounter), trackLeftEye=trackLeftEye, trackRightEye=trackRightEye)
    LTcomment('calibration %d'%(calibrationcounter))
    calibrationcounter += 1

    livetrack_filename = '%sdist_%s_%s_%d_block%d.csv'%(eyetrack_path, hemifield, expInfo['ID'], run_no, blockcounter)

    LiveTrack.SetDataFilename(livetrack_filename)

while any(stairs_ongoing):

    increment = True

    # ## choose staircase
    which_stair = random.choice(list(compress([0,1,2,3,4,5,6,7], stairs_ongoing)))

    # ## set trial
    if position[which_stair] == []:
        random.shuffle(pos_array)
        position[which_stair] = pos_array[:]
    pos = position[which_stair].pop()

    # get position of target and foil:
    # pos         = pos_array[pos]
    target_pos  = positions[pos[0]]
    foil_pos    = positions[pos[1]]

    # shift = [random.choice([-.67, -.33, 0, .33, .67]),random.choice([-.5, -.25, 0, .25, .5])]
    shift = random.sample([-1, -.5, 0, .5, 1],2)
    dif = intervals[cur_int[which_stair]] * foil_type[which_stair]
    which_first = random.choice(['Targ', 'Foil'])


    if which_first == 'Targ':
        point_1.pos = pol2cart(target_pos[0][0], target_pos[0][1]       + shift[0])
        point_2.pos = pol2cart(target_pos[1][0], target_pos[1][1]       + shift[0])
        point_3.pos = pol2cart(foil_pos[0][0],   foil_pos[0][1]         + shift[1])
        point_4.pos = pol2cart(foil_pos[1][0],   foil_pos[1][1]   + dif + shift[1])

        if eye[which_stair] == 'left':
            point_1.fillColor = col_left
            point_2.fillColor = col_left
            point_3.fillColor = col_left
            point_4.fillColor = col_left
        else:
            point_1.fillColor = col_right
            point_2.fillColor = col_right
            point_3.fillColor = col_right
            point_4.fillColor = col_right

    else:
        point_3.pos = pol2cart(target_pos[0][0], target_pos[0][1]       + shift[0])
        point_4.pos = pol2cart(target_pos[1][0], target_pos[1][1]       + shift[0])
        point_1.pos = pol2cart(foil_pos[0][0],   foil_pos[0][1]         + shift[1])
        point_2.pos = pol2cart(foil_pos[1][0],   foil_pos[1][1]   + dif + shift[1])


        if eye[which_stair] == 'left':
            point_3.fillColor = col_left
            point_4.fillColor = col_left
            point_1.fillColor = col_left
            point_2.fillColor = col_left
        else:
            point_3.fillColor = col_right
            point_4.fillColor = col_right
            point_1.fillColor = col_right
            point_2.fillColor = col_right

    if any([trackLeftEye, trackRightEye]):
        LiveTrack.StartTracking()

        # LiveTrack.SetDataComment('start trial %d'%(trialcounter))
        LTcomment('start trial %d'%(trialcounter))

        # LiveTrack.SetDataComment('checking drift')
        LTcomment('checking drift')

        # get gaze drift:
        drift_check = LTrefix(cfg, fixpoint=[0,0])
        gaze_drifts = drift_check['offsets']
        print(gaze_drifts)

        if drift_check['calibrated']:
            LTstoreCalibration('%scalibration_%d.json'%(eyetrack_path, calibrationcounter), trackLeftEye=True, trackRightEye=True)
            LTcomment('calibration %d'%(calibrationcounter))
            calibrationcounter += 1

        # LiveTrack.SetDataComment('drift: %0.4f, %0.4f, %0.4f, %0.4f'%(gaze_drifts[0][0], gaze_drifts[0][1], gaze_drifts[1][0], gaze_drifts[1][1]))
        if trackLeftEye:
            LTcomment('left drift: %0.4f, %0.4f'%(gaze_drifts[0][0], gaze_drifts[0][1]))
        if trackRightEye:
            LTcomment('right drift: %0.4f, %0.4f'%(gaze_drifts[0+trackLeftEye][0], gaze_drifts[0+trackLeftEye][1]))

    ## trial
    trial_clock.reset()

    waiting_for_response = True
    tick_played = False

    comment_pair1_on = False
    comment_pair2_on = False
    comment_pair1_off = False
    comment_pair2_off = False

    comment_response_on = False

    hiFusion.resetProperties()
    loFusion.resetProperties()

    if any([trackLeftEye, trackRightEye]):
        # LiveTrack.SetDataComment('showing fixation cross')
        LTcomment('fixation')
    
        calibration_triggered = False

    while waiting_for_response:

        blindspot.draw()
        hiFusion.draw()
        loFusion.draw()

        if any([trackLeftEye, trackRightEye]):

            # check if gaze is close to drifted fixation:
            if trial_clock.getTime() > .2 and trial_clock.getTime() < 1.4: # was 1.8
                data = LiveTrack.GetLastResult()
                if trackLeftEye and data.Tracked:
                    lep  = [ data.GazeX-fixation.pos[0],      data.GazeY-fixation.pos[1]      ]
                    
                    l_dist = ((lep[0]-gaze_drifts[0][0])**2 + (lep[1]-gaze_drifts[0][1])**2)**0.5
                    
                    if (l_dist > 1.5):
                        print(gaze_drifts)
                        print('left eye:')
                        print(lep)
                        print(l_dist)
                        LTcomment('fixation broken')
                        position[which_stair] = position[which_stair] + [pos]
                        increment = False
                        resp = 'no_fixation'
                        correct = 'no_fixation'
                        reversal = 'no_fixation'
                        waiting_for_response = False

                if trackRightEye and data.TrackedRight:
                    rep  = [ data.GazeXRight-fixation.pos[0], data.GazeYRight-fixation.pos[1] ]

                    r_dist = ((rep[0]-gaze_drifts[0+trackLeftEye][0])**2 + (rep[1]-gaze_drifts[0+trackLeftEye][1])**2)**0.5

                    if (r_dist > 1.5):
                        print(gaze_drifts)
                        print('right eye:')
                        print(rep)
                        print(r_dist)
                        LTcomment('fixation broken')
                        position[which_stair] = position[which_stair] + [pos]
                        increment = False
                        resp = 'no_fixation'
                        correct = 'no_fixation'
                        reversal = 'no_fixation'
                        waiting_for_response = False

                    # eye_tracked.append(True)
                    # eye_positions.append([lep, rep])
                    # eye_distances.append(np.mean([np.sqrt(lep[0]**2 + lep[1]**2),np.sqrt(rep[0]**2 + rep[1]**2)]))
                
                else:
                    # eye_tracked.append(False)
                    pass



        if trial_clock.getTime() < .1:
            fixation.lineColor = 'white'
        else:
            fixation.lineColor = 'black'

        fixation.draw()

        if any([trackLeftEye, trackRightEye]):
            if trial_clock.getTime() > .1:
                keys = event.getKeys(keyList=['q','Q'])
                if len(keys):
                    if 'q' in keys or 'Q' in keys:
                        position[which_stair] = position[which_stair] + [pos]
                        increment = False
                        resp = 'recalibrate'
                        correct = 'recalibrate'
                        reversal = 'recalibrate'
                        waiting_for_response = False
                        LTcomment('calibration requested')
                        calibration_triggered = True


        if trial_clock.getTime() > .2 and trial_clock.getTime() < 1.0: # was 1.2
            point_1.draw()
            point_2.draw()
            if comment_pair1_on == False:
                # LiveTrack.SetDataComment('pair 1 on')
                if any([trackLeftEye, trackRightEye]):
                    LTcomment('pair 1 on')
                comment_pair1_on = True

        if trial_clock.getTime() >= 1.0:   # was 1.2
            if comment_pair1_off == False:
                # LiveTrack.SetDataComment('pair 1 off')
                if any([trackLeftEye, trackRightEye]):
                    LTcomment('pair 1 off')
                comment_pair1_off = True 

        if trial_clock.getTime() > .6 and trial_clock.getTime() < 1.4: # was 0.8 and 1.8
            point_3.draw()
            point_4.draw()
            if comment_pair2_on == False:
                # LiveTrack.SetDataComment('pair 2 on')
                if any([trackLeftEye, trackRightEye]):
                    LTcomment('pair 2 on')
                comment_pair2_on = True

        if trial_clock.getTime() >= 1.4:  # was 1.8
            if comment_pair2_off == False:
                # LiveTrack.SetDataComment('pair 2 off')
                if any([trackLeftEye, trackRightEye]):
                    LTcomment('pair 2 off')
                comment_pair2_off = True 

        
        # if trial_clock.getTime() > .8 and not(tick_played):
        #     tick.play()
        #     tick_played = True

        if trial_clock.getTime() > 1.4: # was 1.2

            if comment_response_on == False:
                # LiveTrack.SetDataComment('response on')
                if any([trackLeftEye, trackRightEye]):
                    LTcomment('waiting for response')
                comment_response_on = True

            fixation.ori = 45
            keys = event.getKeys(keyList=['space','escape','left','right'])
            if len(keys):
                
                if 'escape' in keys:
                    respFile.write("Run manually ended at " + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "!")
                    abort = True
                    waiting_for_response = False
                    if any([trackLeftEye, trackRightEye]):
                        LTcomment('run aborted')
                if 'space' in keys:                
                    position[which_stair] = position[which_stair] + [pos]
                    increment = False
                    resp = 'abort'
                    correct = 'abort'
                    reversal = 'abort'
                    waiting_for_response = False
                    if any([trackLeftEye, trackRightEye]):
                        LTcomment('trial aborted')
                if 'left' in keys:
                    resp = 1
                    reversal = 0
                    correct = ((which_first == 'Targ') == (dif > 0)) == (True)
                    waiting_for_response = False
                if 'right' in keys:
                    resp = 2
                    reversal = 0
                    correct = ((which_first == 'Targ') == (dif > 0)) == (False)
                    waiting_for_response = False
            
            if waiting_for_response == False:
                # prevent stuck-key syndrome responding to 100 trials:
                event.clearEvents(eventType='keyboard')
                # LiveTrack.SetDataComment('response given')
                if any([trackLeftEye, trackRightEye]):
                    LTcomment('response given')

        win.flip()
    
    # LiveTrack.SetDataComment('end of trial')
    if any([trackLeftEye, trackRightEye]):
        LTcomment('end of trial')
        LiveTrack.StopTracking()

    if abort:
        break

    # ## response
    # fixation.draw()
    # win.flip()
    # k = ['wait']
    # while k[0] not in ['escape', 'space', 'left', 'right']:
    #     k = event.waitKeys()
    # if k[0] in ['escape']:
    #     respFile.write("Run manually ended at " + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "!")
    #     break
    # elif k[0] in ['space']:
    #     # position[which_stair] = position[which_stair] + [pos]
    #     # increment = False
    #     stimuli += [stimulus]
    #     resp = 'abort'
    #     correct = 'abort'
    #     reversal = 'abort'
    # else:
    #     resp = 1 if k[0] == 'left' else 2
    #     reversal = 0
    #     correct = ((which_first == 'Targ') == (dif > 0)) == (k[0] == 'left')

    if increment:

        # store which pair was chosen for this staircase:
        pair_chosen = [['Foil','Targ'],['Targ','Foil']][which_first == 'Targ'][resp-1]
        pairs_chosen[which_stair].append(pair_chosen)

        # this could be changed to a different number, but we want to be done fast!
        n_stair_trials = 1

        # see if we can check if a step needs to be taken:
        if len(pairs_chosen[which_stair]) >= n_stair_trials:
            # there are at least N entries
            lastfew = pairs_chosen[which_stair][-n_stair_trials:]
            # check if we need to take a step:
            if len(set(lastfew)) == 1:
                # all responses the same: move toward PSE
                step_direction = {'Targ':-1, 'Foil':1}[lastfew[0]]
                step_direction = step_direction * (-1 * foil_type[which_stair])
                
                # if people keep indicating the Target is shorter,
                # we need to make the Target longer, which is at the start of the interval list
                # and, vice versa, if people keep indicating the Foil

                # if the step is in a different direction from before, we log it as a reversal:
                if step_direction != last_step_directions[which_stair]:
                    revs[which_stair] += 1
                # and we record the info 
                last_step_directions[which_stair] = step_direction

                # then we make the actual step:
                if step_direction == 1:
                    cur_int[which_stair] = min(max_int_idx, cur_int[which_stair] + 1)
                if step_direction == -1:
                    cur_int[which_stair] = max(0, cur_int[which_stair] - 1)

                # also empty the list of chosen pairs:

                pairs_chosen[which_stair] = []
                
                # because this needs to be decided within steps, not between.. right?

        # these lines from the previous implementation, should still be used:
        trial[which_stair] = trial[which_stair] + 1

        # changing this so that either criteria ends the staircase:
        # stairs_ongoing[which_stair] = revs[which_stair] >= nRevs and trial[which_stair] < nTrials

        # requires both criteria before stopping a staircase:
        stairs_ongoing[which_stair] = revs[which_stair] < nRevs or trial[which_stair] < nTrials

        # correct = ((which_first == 'Targ') == (dif > 0)) == (k[0] == 'left')

        # ## update staircase (which direction, is there a reversal?)
        # reversal = False
        # resps[which_stair] = resps[which_stair] + [correct]
        # if resps[which_stair][-2:] == [True, True] and abs(dif) > 1.5:
        #     direction[which_stair] = 1 if intervals[cur_int[which_stair]] > 0 else -1
        # elif resps[which_stair][-2:] == [False, True]:
        #     reversal = True
        #     revs[which_stair] = revs[which_stair] + 1
        #     direction[which_stair] = 1 if intervals[cur_int[which_stair]] > 0 else -1
        # elif resps[which_stair][-2:] == [True, False]:
        #     reversal = True
        #     revs[which_stair] = revs[which_stair] + 1
        #     direction[which_stair] = -1 if intervals[cur_int[which_stair]] > 0 else 1

        # ## increment/update
        # cur_int[which_stair] = max(min(cur_int[which_stair] + direction[which_stair], len(intervals) - 1), 0)
        # trial[which_stair] = trial[which_stair] + 1
        # stairs_ongoing[which_stair] = revs[which_stair] <= nRevs or trial[which_stair] < nTrials

    fixation.ori = 0

    ## print trial
    print(blockcounter,
          trialcounter,
          resp,
          pos[0],
          pos[1],
          tar,
          dif,
          which_first,
          correct,
          reversal,
          foil_type[which_stair],
          eye[which_stair],
          which_stair)
    respFile.write('\t'.join(map(str, [ blockcounter,
                                        trialcounter,
                                        resp,
                                        pos[0],
                                        pos[1],
                                        tar,
                                        dif,
                                        which_first,
                                        correct,
                                        reversal,
                                        foil_type[which_stair],
                                        eye[which_stair],
                                        which_stair])) + "\n")

    blocktrialcounter += 1
    trialcounter += 1

    if any([trackLeftEye, trackRightEye]):
        if calibration_triggered:

            LiveTrack.StopTracking()

            LTcal(cfg, trackLeftEye=trackLeftEye, trackRightEye=trackRightEye)

            LiveTrack.SetResultsTypeCalibrated()

            LTstoreCalibration('%scalibration_%d.json'%(eyetrack_path, calibrationcounter), trackLeftEye=trackLeftEye, trackRightEye=trackRightEye)
            LTcomment('calibration %d'%(calibrationcounter))
            calibrationcounter += 1

            LiveTrack.StartTracking()




    if (blocktrialcounter >= 40):

        # do a break
        visual.TextStim(win,'time for a break!\n\npress SPACE to end your break\n(and do a new calibration)', height = letter_height, color = 'white').draw()
        win.flip()


        # LiveTrack.SetDataComment('end of block')
        if any([trackLeftEye, trackRightEye]):
            LTcomment('end of block')
            LiveTrack.CloseDataFile()


        waiting_for_break_end = True

        while waiting_for_break_end:
            keys = event.getKeys(keyList=['space','escape'])
            if len(keys):
                if 'space' in keys:
                    waiting_for_break_end = False

        blocktrialcounter = 0

        # redo calibration:
        if any([trackLeftEye, trackRightEye]):
            LiveTrack.StopTracking()

            LTcal(cfg, trackLeftEye=trackLeftEye, trackRightEye=trackRightEye)

            LiveTrack.SetResultsTypeCalibrated()

            LTstoreCalibration('%scalibration_%d.json'%(eyetrack_path, calibrationcounter), trackLeftEye=trackLeftEye, trackRightEye=trackRightEye)
            LTcomment('calibration %d'%(calibrationcounter))
            calibrationcounter += 1

            livetrack_filename = '%sdist_%s_%s_%d_block%d.csv'%(eyetrack_path, hemifield, expInfo['ID'], run_no, blockcounter)
            LiveTrack.SetDataFilename(livetrack_filename)

            LiveTrack.StartTracking()

        blockcounter += 1

respFile.close()

## last screen
blindspot.autoDraw = False
visual.TextStim(win,'Run ended.\n ', height = letter_height, color = 'black').draw()
win.flip()

if any([trackLeftEye, trackRightEye]):
    LTcomment('end of run')
    LiveTrack.CloseDataFile()

    # close LiveTrack:
    LiveTrack.Close()

print(datetime.datetime.now().strftime("%Y-%m-%d-%H-%M"))

visual.TextStim(win,'Run ended.\npress SPACE or ESCAPE to exit', height = letter_height, color = 'black').draw()
k = ['wait']
while k[0] not in ['escape','space']:
    k = event.waitKeys()
if k[0] in ['escape','space']:
    win.close()
    core.quit()



