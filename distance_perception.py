 #!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Distance comparison across blind spot
TWCF IIT vs PP experiment 2a piloting
Authors: Clement Abbatecola, Belén María Montabes de la Cruz
    Code version:
        2.0 # 2024/02/12    Final common version before eye tracking
"""

import sys, os
sys.path.append(os.path.join('..', 'EyeTracking'))
from EyeTracking import localizeSetup, EyeTracker

import psychopy
from psychopy import core, visual, gui, data, event
from psychopy.tools.coordinatetools import pol2cart, cart2pol
import numpy as np
import random, datetime
from glob import glob
from itertools import compress
# from fusion_stim import fusionStim

######
#### Initialize experiment
######

def doDistanceTask(ID=None, side=None):

    ## parameters
    nRevs   = 10   #
    nTrials  = 30  # at least 10 reversals and 30 trials for each staircase (~ 30*8 staircases = 250 trials)
    letter_height = 40

    ## path
    # main_path = "C:/Users/clementa/Nextcloud/project_blindspot/blindspot_eye_tracker"

    if any([ID == None, side == None]):
        
        expInfo = {}

        if ID == None:
            expInfo['ID'] = ''
        if side == None:
            expInfo['side'] = ['RH','LH']

        dlg = gui.DlgFromDict(expInfo, title='Infos', screen=0)
        
        if ID == None:
            ID = expInfo['ID']
        if side == None:
            side = expInfo['side']

    data_path = "../data/distance/" + ID + "/"
    os.makedirs(data_path, exist_ok=True)

    ## open file for behavioral data:
    x = 1
    filename = 'dist_' + side + '_' + ID.lower() + '_'
    while (filename + str(x) + '.txt') in os.listdir(data_path): x += 1
    respFile = open(data_path + filename + str(x) + '.txt','w')

    respFile.write(''.join(map(str, ["Start: \t" + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "\n"])))
    # respFile.write("Resp\tTarg_loc\tFoil_loc\tTarg_len\tDifference\tWhich_first\tTarg_chosen\tReversal\tFoil_type\tEye\tGaze_out\tStair\tTrial\n")
    respFile.write("Resp\tTarg_loc\tFoil_loc\tTarg_len\tDifference\tWhich_first\tTarg_chosen\tReversal\tEye\tGaze_out\tStair\tTrial\n")
    print(datetime.datetime.now().strftime("%Y-%m-%d-%H-%M"))
    # print("Resp", "Targ_loc", "Foil_loc", "Targ_len", "Difference", "Which_first", "Targ_chosen", "Reversal", "Foil_type", "Eye", "Gaze_out", "Stair")
    print("Resp", "Targ_loc", "Foil_loc", "Targ_len", "Difference", "Which_first", "Targ_chosen", "Reversal", "Eye", "Gaze_out", "Stair")

    ## blindspot parameters
    bs_file = open(glob("../data/mapping/" + ID + "_" + side + "_blindspot*.txt")[-1],'r')
    bs_param = bs_file.read().replace('\t','\n').split('\n')
    bs_file.close()

    # spot is not sided
    spot_cart = eval(bs_param[1])
    spot = cart2pol(spot_cart[0], spot_cart[1])
    spot_size = eval(bs_param[3])

    # scale blindspot marker
    # we want half the surface? maybe a little more?
    bs_scale = 0.75
    bsf = bs_scale**0.5
    marker_size = [spot_size[0] * bsf, spot_size[1] * bsf]
    
    '''
    distance of reference between dots (target)
    => width of blindspot + 2 (dot width, padding) + 2 (to account for a max jitter of 1 on either side)
    '''
    tar =  spot_size[0] + 2 + 2

    # size of blind spot + 2 (dot width, padding)
    if spot_cart[1] < 0:
        ang_up = (cart2pol(spot_cart[0], spot_cart[1] - spot_size[1])[0] - spot[0]) + 2
    else:
        ang_up = (cart2pol(spot_cart[0], spot_cart[1] + spot_size[1])[0] - spot[0]) + 2

    ## colour (eye) parameters
    col_file = open(glob("../data/color/" + ID + "_col_cal*.txt")[-1],'r')
    col_param = col_file.read().replace('\t','\n').split('\n')
    col_file.close()

    col_back  = eval(col_param[1])
    col_left  = eval(col_param[3])
    col_right = eval(col_param[5])
    col_both  = [col_left[0], -1, col_right[2]] # do we even need this?


    # this could be smarter (IP address based?) but it doesn't matter:
    if os.sys.platform == 'linux':
        location = 'toronto'
    else:
        location = 'glasgow'

    glasses = 'RG'
    trackEyes = [True, True]

    eyetrackfolder =  "../data/distance/" + ID + "/eyetracking/"
    os.makedirs(eyetrackfolder, exist_ok=True)

    setup = localizeSetup(location=location, glasses=glasses, trackEyes=trackEyes, filefolder=eyetrackfolder) # eye-tracking data stored in sub-folder

    cfg = {}
    cfg['hw'] = setup

    pyg_keyboard = key.KeyStateHandler()
    cfg['hw']['win'].winHandle.push_handlers(pyg_keyboard)


    fixation = visual.ShapeStim(cfg['hw']['win'], vertices = ((0, -2), (0, 2), (0,0), (-2, 0), (2, 0)), lineWidth = 4, units = 'pix', size = (10, 10), closeShape = False, lineColor = 'black')

    ## instructions
    visual.TextStim(cfg['hw']['win'],'Troughout the experiment you will fixate at a white cross that will be located at the right hand side of the screen.   \
    It is important that you fixate on this cross at all times.\n\n You will be presented with pairs of dots. You will have to indicate which dots were closer together.\n\n Left arrow = first pair of dots were closer together.\
    \n\n Right arrow = second pair of dots were closer together.\n\n\n Press the space bar to start the experiment.', height = letter_height,wrapWidth=1200, color = 'black').draw()
    cfg['hw']['win'].flip()
    k = ['wait']
    while k[0] not in ['q','space']:
        k = event.waitKeys()
    if k[0] in ['q']:
        win.close()
        core.quit()


    ######
    #### Prepare stimulation
    ######

    ## stimuli
    # point_1 = visual.Circle(cfg['hw']['win'], radius = .5, pos = pol2cart(00, 3), units = 'deg', fillColor = col_both, lineColor = None)
    # point_2 = visual.Circle(cfg['hw']['win'], radius = .5, pos = pol2cart(00, 6), units = 'deg', fillColor = col_both, lineColor = None)
    # point_3 = visual.Circle(cfg['hw']['win'], radius = .5, pos = pol2cart(45, 3), units = 'deg', fillColor = col_both, lineColor = None)
    # point_4 = visual.Circle(cfg['hw']['win'], radius = .5, pos = pol2cart(45, 6), units = 'deg', fillColor = col_both, lineColor = None)
    # TMI, probably don't even need lineColor (default value: False)
    point_1 = visual.Circle(cfg['hw']['win'], lineColor = None)
    point_2 = visual.Circle(cfg['hw']['win'], lineColor = None)
    point_3 = visual.Circle(cfg['hw']['win'], lineColor = None)
    point_4 = visual.Circle(cfg['hw']['win'], lineColor = None)


    blindspot = visual.Circle( win       = cfg['hw']['win'], 
                               pos       = spot_cart, 
                               fillColor = {'LH':col_left, 'RH':col_right}[side], 
                               lineColor = None, 
                               size      = marker_size)
    # blindspot.autoDraw = True 

    # blindspot marker should not be on autodraw (no stimulus should ever be on autodraw, that's a recipe for mistakes)


    ## prepare trials
    if side == 'LH':
        positions = {
            "top": [(spot[0] - ang_up, spot[1] - tar/2), (spot[0] - ang_up, spot[1] + tar/2)],
            "mid": [(spot[0]         , spot[1] - tar/2), (spot[0]         , spot[1] + tar/2)],
            "bot": [(spot[0] + ang_up, spot[1] - tar/2), (spot[0] + ang_up, spot[1] + tar/2)],
        }

    if side == 'RH':
        positions = {
            "top": [(spot[0] + ang_up, spot[1] - tar/2), (spot[0] + ang_up, spot[1] + tar/2)],
            "mid": [(spot[0]         , spot[1] - tar/2), (spot[0]         , spot[1] + tar/2)],
            "bot": [(spot[0] - ang_up, spot[1] - tar/2), (spot[0] - ang_up, spot[1] + tar/2)],
        }

    # First column is target, second column is foil, those are the valid trials in the proportion we want
    pos_array = [["mid", "top"],
                 ["mid", "bot"],
                 ["top", "bot"],
                 ["bot", "top"]]

    # what is this for?
    pos_array_bsa = pos_array[0:2]
    pos_array_out = pos_array[2:4]

    # random.shuffle(pos_array)

    ######
    #### Prepare eye tracking
    ######

    # initialize eye-tracker + gaze ok region etc.
    cfg['hw']['tracker'].initialize()
    cfg['hw']['tracker'].calibrate()
    cfg['hw']['tracker'].startcollecting()


    fixation.draw()
    cfg['hw']['win'].flip()

    k = event.waitKeys()
    if k[0] in ['q']:
        respFile.close()
        cfg['hw']['tracker'].shutdown()
        cfg['hw']['win'].close()
        core.quit()

    ######
    #### Staircase
    ######

    trial_clock = core.Clock()

    intervals = [3.5, 3, 2.5, 2, 1.5, 1, .5, 0, -.5, -1, -1.5, -2, -2.5, -3, -3.5]
    max_int_idx = len(intervals) - 1

    # foil_type  = [1, -1]                                      * 4
    cur_int    = [0, max_int_idx]                             * 4 # we can do without the foil_type parameter if cur_int can be the minimum (0) OR the maximum...
    eye        = ['left', 'left', 'right', 'right']           * 2
    pos_arrays = [pos_array_bsa[:]] * 4 + [pos_array_out[:]]  * 4 

    position       = [[]]    * 8
    trial_stair    = [0]     * 8
    revs           = [0]     * 8
    direction      = [1]     * 8
    # cur_int        = [0]     * 8 # weren't we supposed to start some at the -3.5 interval? now they all start at 3.5... foil_type gets multiplied by this?
    reversal       = [False] * 8
    resps          = [[]]    * 8
    stairs_ongoing = [True]  * 8

    trial = 1
    abort = False
    recalibrate = False

    while any(stairs_ongoing):

        # do a trial
        fixating = False
        while not fixating:
            if cfg['hw']['tracker'].gazeInFixationWindow():
                fixating = True
            else:
                if cfg['hw']['tracker'].calibrate():
        
        # good to go:
        increment = True

        ## choose staircase
        which_stair = random.choice(list(compress([x for x in range(len(stairs_ongoing))], stairs_ongoing)))

        ## set trial
        if position[which_stair] == []:
            random.shuffle(pos_arrays[which_stair])
            position[which_stair] = pos_arrays[which_stair][:]
        pos = position[which_stair].pop()

        shift = random.sample([-1, -.5, 0, .5, .1], 2)
        dif = intervals[cur_int[which_stair]] #* foil_type[which_stair]
        which_first = random.choice(['Targ', 'Foil'])

        target_pos  = positions[pos[0]]
        foil_pos    = positions[pos[1]]

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

        
        cfg['hw']['fusion']['hi'].resetProperties()
        cfg['hw']['fusion']['lo'].resetProperties()

        # record trial number:
        cfg['hw']['tracker'].comment('start trial %d'%(trial))

        # pre-trial fixation
        cfg['hw']['tracker'].waitForFixation()

        # new random patterns:
        cfg['hw']['fusion']['hi'].resetProperties()
        cfg['hw']['fusion']['lo'].resetProperties()

        # now we start the clock:
        trial_clock.reset()

        response     = None
        recalibrate  = False

        while response == None:

            t = trial_clock.getTime()

            if 0.1 <= trial_clock.getTime() < 0.9:
                point_1.draw()
                point_2.draw()
            if 0.5 <= trial_clock.getTime() < 1.3:
                blindspot.draw()
                point_3.draw()
                point_4.draw()
            
            if trial_clock.getTime() < 1.3:
                blindspot.draw()
                k = event.getKeys(['q','r','space'])
                if not cfg['hw']['tracker'].gazeInFixationWindow():
                    response    = 'abort'
                    # recalibrate = True   # maybe only if this happens a few times in a row? or triggered by experimenter...
            else:
                k = event.getKeys(['q','r','space','left','right'])
                fixation.ori += 45
            
            cfg['hw']['fusion']['hi'].draw()
            cfg['hw']['fusion']['lo'].draw()
            fixation.draw()

            cfg['hw']['win'].flip()

            if k and 'q' in k:
                respFile.write("Run manually ended at " + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "!")
                return(0)

            if k and 'space' in k:
                response = 'abort'

            if k and 'r' in k:
                response = 'abort'
                recalibrate = True

            if k and k[0] in ['left','right']:
                response = 1 if k[0] == 'left' else 2
            
                
        fixation.ori -= 45
        
        if abort:
            position[which_stair] = position[which_stair] + [pos] # what does this do?
            increment = False
            resp        = 'abort'
            targ_chosen = 'abort'
            reversal    = 'abort'
            
        # else:
        
        #     ## dealing with auto-aborted trials
        
        #     # auto recalibrate if no initial fixation
        #     if recalibrate:
        #         # recalibrate = False
        #         # visual.TextStim(win,'Calibration...', color = 'black', units = 'deg', pos = (0,-2)).draw()
        #         # fixation.draw()
        #         # win.flip()
        #         # k = event.waitKeys()
        #         # if k[0] in ['q']:
        #         #     respFile.write("Run manually ended at " + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "!")
        #         #     break

        #         # stoprecording()?
        #         cfg['hw']['tracker'].calibrate()
        #         # startrecording()?
                    
        #         #!!# calibrate
                
        #         fixation.draw()
        #         cfg['hw']['win'].flip()
        #         k = event.waitKeys()
        #         if k[0] in ['q']:
        #             respFile.write("Run manually ended at " + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "!")
        #             break
            
            # changing fixation to signify gaze out, restart with 'up' possibily of break and manual recalibration 'r' 
            # this just makes the experiment longer and more frustrating?
            # else:
            #     cfg['hw']['fusion']['hi'].draw()
            #     cfg['hw']['fusion']['lo'].draw()
            #     visual.TextStim(cfg['hw']['win'], '#', height = letter_height, color = 'black').draw()
            #     cfg['hw']['win'].flip()
            #     k = ['wait']
            #     while k[0] not in ['q', 'up', 'r']:
            #         k = event.waitKeys()
                # if k[0] in ['q']:
                #     respFile.write("Run manually ended at " + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "!")
                #     break
        
                # manual recalibrate
                if k[0] in ['r']:
                    visual.TextStim(cfg['hw']['win'],'Calibration...', color = 'black', units = 'deg', pos = (0,-2)).draw()
                    fixation.draw()
                    cfg['hw']['win'].flip()
                    k = event.waitKeys()
                    if k[0] in ['q']:
                        respFile.write("Run manually ended at " + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "!")
                        break

                    #!!# calibrate

                    fixation.draw()
                    cfg['hw']['win'].flip()
                    k = event.waitKeys()
                    if k[0] in ['q']:
                        respFile.write("Run manually ended at " + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M") + "!")
                        break
                
            position[which_stair] = position[which_stair] + [pos]
            increment = False
            resp = 'auto_abort'
            targ_chosen = 'auto_abort'
            reversal = 'auto_abort'
        
        if increment:
            '''
            which_first == 'Targ'          => was target first? (True/False)
            dif > 0                        => was target smaller? (True/False)
            k[0] == 'left'                 => was first chosen? (True/False)
            target first == target smaller => was first smaller? (True/False)
            first smaller == first chosen  => was smaller chosen? (True/False)
            
            (which_first == 'Targ') == (k[0] == 'left') => was target chosen?
            '''
            
            targ_chosen = (which_first == 'Targ') == (k[0] == 'left')

            ## update staircase (which direction, is there a reversal?)
            reversal = False
            resps[which_stair] = resps[which_stair] + [targ_chosen]
            if (len(resps[which_stair]) > 1) and (resps[which_stair][-2] != resps[which_stair][-1]):
                reversal = True
                revs[which_stair] = revs[which_stair] + 1
                direction[which_stair] *= -1
            
            ## increment/update
            cur_int[which_stair] = max(min(cur_int[which_stair] + direction[which_stair], len(intervals) - 1), 0)
            trial_stair[which_stair] = trial_stair[which_stair] + 1
            stairs_ongoing[which_stair] = revs[which_stair] <= nRevs or trial_stair[which_stair] < nTrials

        

        ## print trial
        print(resp,
            pos[0],
            pos[1],
            tar,
            dif,
            which_first,
            targ_chosen,
            reversal,
            # foil_type[which_stair],
            eye[which_stair],
            gaze_out,
            which_stair)
        respFile.write('\t'.join(map(str, [resp,
                                        pos[0],
                                        pos[1],
                                        tar,
                                        dif,
                                        which_first,
                                        targ_chosen,
                                        reversal,
                                        # foil_type[which_stair],
                                        eye[which_stair],
                                        gaze_out,
                                        which_stair,
                                        trial])) + "\n")
        trial += 1

    if not any(stairs_ongoing):
        print('run ended properly!')

    respFile.close()
    print(datetime.datetime.now().strftime("%Y-%m-%d-%H-%M"))
    blindspot.autoDraw = False

    #!!# close eye-tracker

    ## last screen
    visual.TextStim(cfg['hw']['win'],'Run ended.', height = letter_height, color = 'black').draw()
    cfg['hw']['win'].flip()
    k = event.waitKeys()

    # need to cleanly shutdown the systems:
    cfg['hw']['tracker'].shutdown()
    cfg['hw']['win'].close()
    core.quit()

