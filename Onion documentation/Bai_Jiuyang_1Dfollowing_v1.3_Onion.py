# Jiuyang Bai 1D following v1.2 Carrot4 (batch output)  
# This version saves data to a temporary string and writes the string at the end of each trial
# This solves the frame rate fluctuation problem
# Read input from one folder and write output to one folder (all subjects in one folder)
# save position and orientation data in one file

# 01/08/2018

import math
import lights # turns on the lights
import random
import csv
#from os.path import exists
import os

# Vizard Imports
import viz
import vizact
import viztracker


# The following libraries are won't work outside of the VENLab without their 
# associated dependencies, but are requried for experiments within the VENLab.
# When experiment is ready to be moved to the VENLab, we'll re-include them.

import oculus
import emergencyWalls

#####################################################################################
# Constants


IPD = viz.input('Please enter IPD:')

# Set to True when ready to have the experiment write data
DATA_COLLECT = True
 
# If program should run practice trials
DO_PRACTICE = True

# If run free walk trials
DO_FREEWALK = True

# If program crashes, start trials here
START_ON_TRIAL = 1
# Total number of trials in experiment
FREEWALK_TRIALS = 4 # 4 trials each session
FREEWALK_SESSIONS = 2 # 1 session before practice 1 session after experiment
PRACTICE_TRIALS = 4
TOTAL_TRIALS = 90    # 3(d0) *  3(dv) * 10(reps) = 90
TRIAL_LENGTH = 12 # 12 seconds


# Used for file naming, currently placeholders
EXPERIMENT = '1Dfollowing_v1.3_Onion'
NICKNAME = 'Onion'
EXPERIMENTER_NAME = 'jiuyangbai'

# Set output directory for writing data
OUTPUT_DIR = '/'.join(['Data', EXPERIMENTER_NAME, EXPERIMENT,(NICKNAME + '_Output')]) + '/'
INPUT_DIR = '/'.join(['Data', EXPERIMENTER_NAME, EXPERIMENT,(NICKNAME + '_Input')]) + '/'
SOUND_DIR = '/'.join(['sounds', EXPERIMENTER_NAME, EXPERIMENT]) + '/'
MODEL_DIR = 'Models/'

# inputFile = 'Data/jiuyangbai/1Dfollowing_v1.1_Carrot2/Subject_00/Input/trial001.csv'


# Orientation constants
POLE_TRIGGER_RADIUS = 0.3 # How close participant must be to home pole
THRESHOLD_THETA = 15 # Maximum angle participant can deviate when looking at orienting pole
ORIENT_TIME = 3 # How long participant must orient onto pole

# The dimension of the room space used for experiment
DIMENSION_X = 9.0 # the length of the shorter side in meter
DIMENSION_Z = 11.0 # the lengtth of the longer side in meter
DIAGONAL = (DIMENSION_X**2 + DIMENSION_Z**2)**(1.0/2)# The length of the diagonal line of the experimental space

ROOM_ANGLE = math.atan(DIMENSION_X/DIMENSION_Z) # the anger (in radian) between the diagonal and the shorter edge of the room

# Home and Orient Pole positions (x,z,y)
HOME_POLE = [[DIMENSION_X/2, 0.0, DIMENSION_Z/2], [-DIMENSION_X/2, 0.0, -DIMENSION_Z/2]]
ORI_POLE = [[DIMENSION_X/2 - DIMENSION_X/3, 0.0, DIMENSION_Z/2 - DIMENSION_Z/3], \
			[-DIMENSION_X/2 + DIMENSION_X/3, 0.0, -DIMENSION_Z/2 + DIMENSION_Z/3]]


# Describe the end-trial trigger line (end line) in the intersense coordinate system
# the end line is perpendicular to walking direction
K = -DIMENSION_X/DIMENSION_Z # The slope of the end line
END_DIS = 2.0 # The distance between end line to home pole position of following trial
# The intercept of the end line for two home poles respectively
B = [-(DIAGONAL/2 - END_DIS) / math.cos(ROOM_ANGLE), (DIAGONAL/2 - END_DIS) / math.cos(ROOM_ANGLE)]

#####################################################################################
# Control Options
# Sets the controls / displays
# Oculus Rift Selection unlikely to work outside the VENLab


# Dialog box asking for type of control and subject number
OCULUS = 'Oculus Rift'
MONITOR = 'PC Monitor'
controlOptions = [OCULUS,MONITOR]
controlType = controlOptions[viz.choose('How would you like to explore? ', controlOptions)]

subject = viz.input('Please enter the subject number:','')
subject = str(subject).zfill(2)



# Use keyboard controls
# Controls:
# q - Strafe L		w - Forward		e - Strafe R
# a - Turn L		s - Back		d - Turn R
#
# y - Face Up		r - Fly Up
# h - Face Down		f - Fly Down
if controlType == MONITOR:
	headTrack = viztracker.Keyboard6DOF()
	link = viz.link(headTrack, viz.MainView)
	headTrack.eyeheight(1.6)
	link.setEnabled(True)
	viz.go()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
# Use Oculus Rift
elif controlType == OCULUS:
	# viz.fullscreen.x is the inital x position when run
	viz.ipd(IPD)	

	#viz.vsync(0)
	#viz.setDisplayMode(0,0,0,60)
	viz.go()
	
	# add intersense tracker
	isense = viz.add('intersense.dle')
	ISTracker = isense.addTracker(port=5001,station=0)
	ISTracker.setEnhancement(2) 
	ISTracker.setSensitivity(4)
	ISTracker.setShockSuppression(2)
	ISTracker.setAccelSensitivity(4)	
	# add Oculus tracker
	OVRTracker = oculus.Rift().getSensor()
	
	# add the virtual tracker, link it to the MainView and set an offset to get the eye position
	virtual = viz.addGroup()
	link = viz.link(virtual, viz.MainView)
	link.preTrans([0, -0.055, -0.073]) # [-, -, -] = [left, down, back]
	
	# tracker initialization flag
	init = False
	
	# yaw correction flag
	corrected = False

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Use Oculus Rift
#elif controlType == OCULUS:
#	hmd = oculus.Rift()
#	viz.setOption('viz.fullscreen.x','0')
#	viz.setOption('viz.fullscreen.y','0')
#	vizconnect.go('vrpnOculusMerge.pyc')
#	hmdName = 'oculus'

viz.clip(.001,1000) # Adjust size of the viewing frustum



######################################################################################################
# Experiment Loop Helper Functions
def goToStage(nextTrialStage):
	global trial_stage
	trial_stage = nextTrialStage
	print 'Going to: ' + trial_stage

	
def endLine(x):
	# takes the x value and gives the z value of the end trigger line defined in the intersense coordinate system
	return x*K + B[trial_num%2]
	
def moveTarget(traj, index):
	
	global posIndex, time, targetTraj
	
	if trial_num%2 == 1: 
		# start from (-, 0, -) corner of the room	
		dz = float(traj[index][1])*math.cos(ROOM_ANGLE) - float(traj[index][0])*math.sin(ROOM_ANGLE)
		dx = float(traj[index][1])*math.sin(ROOM_ANGLE) + float(traj[index][0])*math.cos(ROOM_ANGLE)
	else:                
		# start from (+, 0, +) corner of the room
		dz = -float(traj[index][1])*math.cos(ROOM_ANGLE) + float(traj[index][0])*math.sin(ROOM_ANGLE)
		dx = -float(traj[index][1])*math.sin(ROOM_ANGLE) - float(traj[index][0])*math.cos(ROOM_ANGLE)
	
	models['targetPole'].setPosition([HOME_POLE[trial_num%2][0]+dx, 0, HOME_POLE[trial_num%2][2]+dz])
	

def countDown(t):
	global time, isStamped, time_stamp
	isTimeUp = False
	if not isStamped:
		time_stamp = time
		isStamped = True
	if time - time_stamp > t:
		isTimeUp = True
	return isTimeUp
	
def writeCSVFile(fileName, data, time):
	strData = [str(round(t,4)) for t in data+[time]]
	file = open(fileName, 'a')
	file.write(','.join(strData)+'\n')
	file.close()
	
def relativeOrientation(pos1, pos2):
	xrel = round(pos2[0]-pos1[0],4)
	zrel = round(pos2[2]-pos1[2],4)
	theta = 0
	if zrel == 0.0 and xrel > 0:
		theta = math.pi/2
	elif zrel == 0.0:
		theta = math.pi/2*3
	else:
		theta = math.atan(round(xrel,4)/round(zrel,4))
		if zrel < 0:
			theta += math.pi
		if zrel > 0 and xrel < 0:
			theta += math.pi*2
	return theta
	
def facing(lookingObjectPosn, lookedAtObjectPosn, lookingObjectYaw, thresholdTheta):
	"""lookingObjectPosn: position of object that is looking
	lookedAtObjectPosn: position of object looked at
	lookingObjectYaw: yaw of the object that is looking (degrees)
	thresholdTheta: viewing angle must be +/- this amount in order to be considered 'looking at' the object. degrees

	return: bool, whether the looking object is facing the looked-at object
	>>> universals.facing([0,0,0],[1,0,5],0,20)
	True
	>>> universals.facing([3,0,3],[1,0,0],210,20)
	True
	"""
	degRelOrientation = 180.0/math.pi*relativeOrientation(lookingObjectPosn, lookedAtObjectPosn) #radians
	degRelOrientation = (degRelOrientation+180)%360-180
	return math.fabs(degRelOrientation-lookingObjectYaw)<thresholdTheta
	
def distance(x,y,a,b):
	return ((x-a)**2+(y-b)**2)**.5
	
def inRadius(pos, center, radius):
	#This method takes in two poadsitions in [x,y,z] form and then returns true if the distance between them is less than the radius given
	if pos == '' or center == '' or radius == '':
		return False
	return (distance(pos[0],pos[2],center[0],center[2]) <= radius)
	
def inRadiusAndFacing ():
	global POLE_TRIGGER_RADIUS, THRESHOLD_THETA, models
	cur_pos = viz.get(viz.HEAD_POS)
	cur_rot = viz.get(viz.HEAD_ORI)
	return inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS) and facing(cur_pos, models['orientPole'].getPosition(), cur_rot[0], THRESHOLD_THETA)
	
#######################################################################################################
# Experiment


# loads condition data
inputFile = INPUT_DIR + NICKNAME + '_subject' + subject + '_conditions.csv'
with open(inputFile, 'rb') as data:
	rows = csv.reader(data)
	conditions = [[value for key, value in enumerate(row)] for row in rows]
condition = ''

viz.clearcolor(0,0.4,1.0) # blue world
models = {}
models['homePole'] = viz.add(MODEL_DIR + 'bluePole.3ds')
models['orientPole'] = viz.add(MODEL_DIR + 'redPole.3ds')
models['targetPole'] = viz.add(MODEL_DIR + 'greenPole.3ds')
models['ground'] = viz.add(MODEL_DIR + 'ground4.osgb')
# Adjust models size
models['homePole'].setScale([0.6,0.45,0.6]) # the original size = [0.4m 3m 0.4m]
#models['targetPole'].setScale([1.0,0.66,1.0]) # the original size = [0.4m 3m 0.4m]
# the transparency
_alpha = 0.0
# Hide loaded models
models['homePole'].visible(viz.OFF)
models['orientPole'].visible(viz.OFF)
models['targetPole'].visible(viz.OFF)


# Sounds
sounds={}
sounds['1intro_preFreewalk'] = viz.addAudio(SOUND_DIR + '1intro_preFreewalk.mp3')
sounds['2preFreewalkHome'] = viz.addAudio(SOUND_DIR + '2preFreewalkHome.mp3')
sounds['3practice'] = viz.addAudio(SOUND_DIR + '3practice.mp3')
sounds['4experimental'] = viz.addAudio(SOUND_DIR + '4experimental.mp3')
sounds['5postFreewalk'] = viz.addAudio(SOUND_DIR + '5postFreewalk.mp3')
sounds['6end'] = viz.addAudio(SOUND_DIR + '6end.mp3')
sounds['begin'] = viz.addAudio(SOUND_DIR + 'begin.mp3')


# Initial data_collection, regardless of DATA_COLLECT
data_collect = False

# Initializes trial_stage, which is the current step of a given trial
# First part is the name of the specific stage (pretrial)
# Second part is practice (01) or experimental trials (02)
# Third part is the stage's position in the order of all stages (01)
goToStage('pretrial_00_01')

# Initializa setting to start free walk trials
if DO_FREEWALK == True:
	is_freewalk = True
else:
	is_freewalk = False
	goToStage('pretrial_01_01')

# Initially setting to start practice trials
if DO_PRACTICE == True:
	is_practice = True
else:
	is_practice = False
	goToStage('pretrial_02_01')
	
# Starts trial count at 1, unless otherwise specifed
if START_ON_TRIAL > 1:
	trial_num = START_ON_TRIAL
	
else:
	trial_num = 1
freewalk_session = 1
# Time counter set to 0
time = 0

instruction = True
isStamped = False

screenshot = 1
videorecord = 1
data_batch = ''

########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################
########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################
########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################
########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################
########################################################################################################################
# Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop ## Master Loop #
########################################################################################################################

def masterLoop(num):
	# global variables within masterLoop
	global DATA_COLLECT, DO_PRACTICE, is_practice, is_freewalk, data_collect, trial_stage, trial_num, \
	freewalk_session, time, time_stamp, cur_pos, posIndex, conditions, condition, K, B, targetTraj, _alpha, \
	isStamped, controlType,instruction, screenshot, videorecord,\
	old_OVROri, old_camOri, init, corrected, virtual, \
	data_batch #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	# Time elapsed since the last run of masterLoop and then added to the global time
	frame_elapsed = viz.getFrameElapsed()
	time += frame_elapsed
	
	viz.setOption('viz.AVIRecorder.maxWidth','1024')
	viz.setOption('viz.AVIRecorder.maxHeight','720')
	viz.setOption('viz.AVIRecorder.fps','30')
	

	if os.path.isfile(OUTPUT_DIR + 'image'+ str(screenshot) +'.bmp') == True:
		screenshot += 1
	if os.path.isfile(OUTPUT_DIR + 'video'+ str(videorecord) +'.avi') == True:
		videorecord += 1
		
	vizact.onkeydown('p', viz.window.screenCapture, OUTPUT_DIR + 'image'+ str(screenshot) +'.bmp')
	vizact.onkeydown('b', viz.window.startRecording, OUTPUT_DIR + 'video'+ str(videorecord) +'.avi')
	vizact.onkeydown('e', viz.window.stopRecording)
	

	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	if controlType == OCULUS:
		if not init:
			try:
				old_OVROri = OVRTracker.getEuler()
				# initial the camera euler by Interesense (Actually just want yaw. Pitch and roll will be replaced by Oculus later)
				old_camOri = ISTracker.getEuler()
				# replace pitch and roll by Oculus because its more sensitive and doesn't give too much error
				old_camOri[1] = old_OVROri[1]
				old_camOri[2] = old_OVROri[2]
			except NameError:
				print('tracker haven''t been initialized')
			else:
				init = True

		else:
		
			# get the current orientation from trackers
			ISOri = ISTracker.getEuler()
			OVROri = OVRTracker.getEuler()
			
			# apply the change of Oculus orientation to the camera
			d_OVROri = [OVROri[0] - old_OVROri[0],\
						OVROri[1] - old_OVROri[1],\
						OVROri[2] - old_OVROri[2]]
						
			camOri = [old_camOri[0] + d_OVROri[0],\
					  old_camOri[1] + d_OVROri[1],\
					  old_camOri[2] + d_OVROri[2]]


			# Update MainView by intersense position and intersense-corrected Oculus orientation
			virtual.setPosition(ISTracker.getPosition())
			virtual.setEuler(camOri)
			# The orientation at this frame become the old orientation for next frame
			old_camOri = camOri
			old_OVROri = OVROri

	#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	
	# Current position and roation of the participant
	cur_pos = viz.get(viz.HEAD_POS)
	cur_rot = viz.get(viz.HEAD_ORI)
	
	#>> Will only work in VENLab <<
	emergencyWalls.popWalls(cur_pos) # Pops up the Emergency Walls when participant is close to physical room edge.


	
	
	##################
	# Begin Freewalk #
	##################
	
	if is_freewalk == True:
		
		
		# Writes Position and Rotation, but only when DATA_COLLECT set to True
		# and trial_stage has set data_collect to True
		if DATA_COLLECT and data_collect:

			
			target_loc = models['targetPole'].getPosition()
			# Position: Target_x, Target_y, Participant_x, Participant_y, time stamp
			data = [target_loc[0],target_loc[2],cur_pos[0], cur_pos[2], cur_rot[0], cur_rot[1], cur_rot[2]]
			strData = [str(round(t,4)) for t in data+[time]]
			strData = ','.join(strData)+'\n'
			data_batch = data_batch + strData
		
		
		#########
		# 00 01 Freewalk Pretrial: sets up practice trial, establishes pole locations
		if trial_stage == 'pretrial_00_01':

			print '> Start Free Walking Session ' + str(freewalk_session) + ' Trial ' + str(trial_num) + ' ----------------------------------'
		
			# Set position of home pole (where participant stands to start trial)			
			if models['homePole'].getVisible() == False:
				models['homePole'].setPosition(HOME_POLE[trial_num%2])
				models['homePole'].alpha(1.0)
				models['homePole'].visible(viz.ON)

			# Set position of orient pole (where participant faces to start trial)
			if models['orientPole'].getVisible() == False:
				models['orientPole'].setPosition(ORI_POLE[(trial_num+1)%2])
				models['orientPole'].visible(viz.ON)
			
			if trial_num == 1 and freewalk_session == 1 and instruction:
				sounds['1intro_preFreewalk'].play()
				instruction = False
			
			if not(trial_num == 1 and freewalk_session == 1) or countDown(53):
				# Move to next stage
				goToStage('orient_00_02')
				instruction = True
				isStamped = False				
				
		#########
		# 00 02 Orienting to Pole: Give time for participant to orient to the pole
		elif trial_stage == 'orient_00_02':
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			if controlType == OCULUS:
				if not corrected and inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS):
					
					
					# Yaw correction
					old_camOri[0] = ISOri[0]
					corrected = True
					print('&&&&&&   rotation corrected 02 02   &&&&&&&&')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			if inRadiusAndFacing():
				if trial_num == 1 and instruction:
					if freewalk_session == 1:
						sounds['2preFreewalkHome'].play()
					elif freewalk_session == 2:
						sounds['5postFreewalk'].play()
					instruction = False
				if not(trial_num == 1) or (freewalk_session == 1 and countDown(16)) or (freewalk_session == 2 and countDown(14)):	
					# Move to stage 3
					goToStage('orient_00_02_wait')
					instruction = True
					isStamped = False
	
				
		#########
		# wait for orientation
		elif (trial_stage == 'orient_00_02_wait'):			
			if countDown(ORIENT_TIME):
				goToStage('inposition_00_03')
				isStamped = False
			if not inRadiusAndFacing():
				isStamped = False
		
		#########
		# 00 03 Freewalk In Position: proceeds once participant is standing on home and facing orient for three seconds
		elif (trial_stage == 'inposition_00_03'):
			print 'Free walk start'
			# Turn off home pole and orient pole
			models['homePole'].visible(viz.OFF)
			models['orientPole'].visible(viz.OFF)
			models['targetPole'].setScale([1.0,0.66,1.0]) # the original size = [0.4m 3m 0.4m]
			models['targetPole'].setPosition(HOME_POLE[(trial_num+1)%2])
			models['targetPole'].visible(viz.ON)
			
			sounds['begin'].play()
	
			# Start to collect data
			data_collect = True
		
			# initialize batch data output
			data_batch = ''
			time = 0
			
			# Move to Stage 4
			goToStage('target_00_04')
			
		#########
		# 00 04 Freewalk: Participants Moves
		elif (trial_stage == 'target_00_04'):
			corrected = False # +++++++++++++++++++++++++++++++++++++++++++++++++++++
				
			# Detects participant location, moves to Stage 5 (Ends Trial) when participant reache the end line
			if (trial_num%2 == 1) and (cur_pos[2] > endLine(cur_pos[0])) or \
				(trial_num%2 == 0)and (cur_pos[2] < endLine(cur_pos[0])):
					
				goToStage('endtrial_00_05')
		
		#########
		# 00 05 End Freewalk Trial: Close out the trial and reset values for next practice trial or start Experiment
		elif trial_stage == 'endtrial_00_05':
			
			# Clears the target pole
			models['targetPole'].visible(viz.OFF)
			
			# save the data of this trial
			fileName = OUTPUT_DIR + NICKNAME + '_freewalk' + \
			'_subj' + subject + '_s' + str(freewalk_session) + '_trial' + str(trial_num).zfill(3) + '.csv'
			file = open(fileName, 'a')
			file.write(data_batch)
			file.close()



			print 'End Freewalk Trial ' + str(trial_num)
			data_collect = False


			# End Check: When trial_num is greater than FREEWALK_TRIALS, end practice and start experiment block
			if trial_num == FREEWALK_TRIALS:
				print '>> End Freewalk Session<<'
	
				if freewalk_session == 2:
					print '>>> End Experiment <<<'
					goToStage('NULL')
					if instruction:						
						sounds['6end'].play()
						instruction = False
				elif freewalk_session == 1:
					goToStage('pretrial_01_01')
					is_freewalk = False
					trial_num = 1
					freewalk_session += 1
			# Returns to Stage 1, resets clock
			else:
				trial_num += 1
				goToStage('pretrial_00_01')				
				instruction = True
		
			
			
	##################
	# Begin practice #
	##################
	
	elif is_practice == True:
		
	
		#########
		# 01 01 Practice Pretrial: sets up practice trial, establishes pole locations
		if trial_stage == 'pretrial_01_01':

			print '> Start Practice Trial ' + str(trial_num) + ' ----------------------------------'
			
			# Loads input files

			inputFile = INPUT_DIR + NICKNAME + '_practice' + str(trial_num) + '.csv'

			with open(inputFile, 'rb') as data:
				rows = csv.reader(data)
				targetTraj = [[value for key, value in enumerate(row)] for row in rows]
	
			# will be used later as the index of targetTraj to update target pole position
			posIndex = 0
	
			# Move to Stage 2
			goToStage('orient_01_02')
					
			
		#########
		# 01 02 Orienting to Pole: Give time for participant to orient to the pole
		elif (trial_stage == 'orient_01_02'):

			# Set position of home pole (where participant stands to start trial)
			if models['homePole'].getVisible() == False:
				models['homePole'].setPosition(HOME_POLE[trial_num%2])
				
			if _alpha < 1.0:
				models['homePole'].alpha(_alpha)
				models['homePole'].visible(viz.ON)
				_alpha += 0.011
	
			# Set position of orient pole (where participant faces to start trial)
			if models['orientPole'].getVisible() == False:
				models['orientPole'].setPosition(ORI_POLE[(trial_num+1)%2])
				models['orientPole'].visible(viz.ON)
	
				
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					# Yaw correction			
			if controlType == OCULUS:	
				if not corrected and inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS):
					
					old_camOri[0] = ISOri[0]
					corrected = True
					print('&&&&&&   rotation corrected 02 02   &&&&&&&&')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			if inRadiusAndFacing():
				if trial_num == 1 and instruction:
					sounds['3practice'].play()
					instruction = False
				if not(trial_num == 1) or countDown(33):
					isStamped = False
					# Move to stage 3
					goToStage('orient_01_02_wait')
					instruction = True
					corrected = False # +++++++++++++++++++++++++++++++++++++++++++++++++++++
					_alpha = 0.0
					# Current time
					time_stamp = time
					
				
				
		#########
		# wait for orientation
		elif (trial_stage == 'orient_01_02_wait'):
			if countDown(ORIENT_TIME):
				goToStage('inposition_01_03')
				isStamped = False
			if not inRadiusAndFacing():
				isStamped = False
				
				
		#########
		# 01 03 Practice In Position: proceeds once participant is standing on home and facing orient for three seconds
		elif trial_stage == 'inposition_01_03':
			
			print 'Practice Target Appears'
			
			# Turn off home pole and orient pole
			models['homePole'].visible(viz.OFF)
			models['orientPole'].visible(viz.OFF)


			# Initializes target pole at home location, which is the reference point for moveTarget function to update target position)
			models['targetPole'].setPosition(HOME_POLE[trial_num%2]) 
			moveTarget(targetTraj, 0) # update the position once to create initial distance before it appears
			size = float(conditions[trial_num][5])
			ratio = size/0.4
			models['targetPole'].setScale([ratio,1,ratio]) # the original size = [0.4m 3m 0.4m]
			models['targetPole'].visible(viz.ON) # target pole appears
			sounds['begin'].play()
			
			# Move to Stage 4
			goToStage('target_01_04')
			time = 0
			
	
		
		#########
		# 01 04 Moving Pole: Target Moves
		elif (trial_stage == 'target_01_04'):

			# Target Moves
			moveTarget(targetTraj, posIndex)
			posIndex = int(time/TRIAL_LENGTH * len(targetTraj))
			
				
			# Detects participant location, moves to Stage 6 (Ends Trial) when participant reache the end line
			if (trial_num%2 == 1) and (cur_pos[2] > endLine(cur_pos[0])) or \
				(trial_num%2 == 0)and (cur_pos[2] < endLine(cur_pos[0])) or \
				posIndex >= len(targetTraj):
				goToStage('endtrial_01_05')
		
		#########
		# 01 05 End Practice Trial: Close out the trial and reset values for next practice trial or start Experiment
		elif trial_stage == 'endtrial_01_05':
			
			# Clears the target pole
			models['targetPole'].visible(viz.OFF)
			
			print 'End Practice Trial ' + str(trial_num)
			

			# End Check: When trial_num is greater than PRACTICE_TRIALS, end practice and start experiment block
			if trial_num >= PRACTICE_TRIALS:
				print '>> End Practice <<'
				goToStage('pretrial_02_01')
				is_practice = False
				trial_num = 1
			# Returns to Stage 1, resets clock
			else:
				trial_num += 1
				goToStage('pretrial_01_01')
			
	####################
	# Begin Experiment #
	####################
	
	elif is_practice == False:
		

		# Writes Position and Rotation, but only when DATA_COLLECT set to True
		# and trial_stage has set data_collect to True
		if DATA_COLLECT and data_collect:
			
			
			# Location of Target Pole
			target_loc = models['targetPole'].getPosition()
			
			# Position: Target_x, Target_y, Participant_x, Participant_y, Participant yaw, pitch, and roll, time stamp
			data = [target_loc[0],target_loc[2],cur_pos[0], cur_pos[2], cur_rot[0], cur_rot[1], cur_rot[2]]
			strData = [str(round(t,4)) for t in data+[time]]
			strData = ','.join(strData)+'\n'
			data_batch = data_batch + strData

			
			# log IPD
			if trial_num == 1:
				file = open(OUTPUT_DIR + NICKNAME + '_subj' + subject + \
				'_IPD_' + str(IPD) + '.txt', 'a')
				file.close()
			
			
		#########
		# 02 01 Experiment Pretrial: sets up trial, establishes pole locations
		if trial_stage == 'pretrial_02_01':
			
			condition = ', '.join([conditions[trial_num][5], conditions[trial_num][3].zfill(4)])
			
			
			# Print start of trial, trial #, and type of trial [pos][speed][turn]
			print '> Start Trial ' + str(trial_num) + ': ' + condition + ' ----------------------------------'
					
			# Loads input files


			inputFile = INPUT_DIR + NICKNAME + '_subject' + subject + '_trial' + str(trial_num).zfill(3) + '.csv'
					
			with open(inputFile, 'rb') as data:
				rows = csv.reader(data)
				targetTraj = [[value for key, value in enumerate(row)] for row in rows]
	
			# will be used later as the index of targetTraj to update target pole position
			posIndex = 0
			
			
			# Move to Stage 2
			goToStage('orient_02_02')
				
				
		#########
		# 02 02 Orienting to Pole: Give time for participant to orient to the pole
		elif (trial_stage == 'orient_02_02'):
			
			## Set position of home pole (where participant stands to start trial)
			if models['homePole'].getVisible() == False:
				models['homePole'].setPosition(HOME_POLE[trial_num%2])
				
			if _alpha < 1.0:
				models['homePole'].alpha(_alpha)
				models['homePole'].visible(viz.ON)
	
				_alpha += 0.016
			# Set position of orient pole (where participant faces to start trial)
			if models['orientPole'].getVisible() == False:
				models['orientPole'].setPosition(ORI_POLE[(trial_num+1)%2])
				models['orientPole'].visible(viz.ON)
			
			
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++		
			if controlType == OCULUS:
				if not corrected and inRadius(cur_pos, models['homePole'].getPosition(), POLE_TRIGGER_RADIUS):
					
					# Yaw correction
					old_camOri[0] = ISOri[0]
					corrected = True
					print('&&&&&&   rotation corrected 02 02   &&&&&&&&')
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			if inRadiusAndFacing():
				if trial_num == 1 and instruction:
					sounds['4experimental'].play()
					instruction = False
				if not(trial_num == 1) or countDown(22):
					# Move to stage 3
					goToStage('orient_02_02_wait')
					corrected = False # +++++++++++++++++++++++++++++++++++++++++++++++++++++
					_alpha = 0.0
					instruction = True
					isStamped = False
				
	
		#########
		# wait for orientation
		elif (trial_stage == 'orient_02_02_wait'):
			if countDown(ORIENT_TIME):
				goToStage('inposition_02_03')
				isStamped = False
			if not inRadiusAndFacing():
				isStamped = False	
	
		
		#########
		# 02 03 In Position: proceeds once participant is standing on home and facing orient
		elif trial_stage == 'inposition_02_03':

			# Turn off home and orient poles
			models['homePole'].visible(viz.OFF)
			models['orientPole'].visible(viz.OFF)
				
			# Initializes target pole at home location, which is the reference point for moveTarget function to update target position)
			models['targetPole'].setPosition(HOME_POLE[trial_num%2]) 
			moveTarget(targetTraj, 0) # update the position once to create initial distance before it appears
			size = float(conditions[trial_num][5])
			ratio = size/0.4
			models['targetPole'].setScale([ratio,1,ratio]) # the original size = [0.4m 3m 0.4m]
			models['targetPole'].visible(viz.ON) # target pole appears
			sounds['begin'].play()
			
			# Turn on data collection for this trial
			data_collect = True
				
			print 'Target Appears'
				
			# Move to Stage 5
			goToStage('target_02_04')
			# initialize batch data output
			data_batch = ''
			time = 0

			
		#########
		# 02 04 Moving Pole: Target Moves
		elif (trial_stage == 'target_02_04'):
			
			# Target Moves			
			moveTarget(targetTraj, posIndex)
			posIndex = int(time/TRIAL_LENGTH * len(targetTraj))
			
			# Detects participant location, moves to Stage 6 (Ends Trial) when participant reache the end line
			if (trial_num%2 == 1) and (cur_pos[2] > endLine(cur_pos[0])) or \
				(trial_num%2 == 0)and (cur_pos[2] < endLine(cur_pos[0])) or \
				posIndex >= len(targetTraj):
				goToStage('endtrial_02_05')
				
			#########
		# 02 05 End Trial: Close out the trial and reset values for next trial
		elif trial_stage == 'endtrial_02_05':
			
			# Clears the target pole
			models['targetPole'].visible(viz.OFF)
			
			# End data collection for this trial
			data_collect = False
			
			# save the data of this trial
			
			fileName = OUTPUT_DIR + NICKNAME + '_subj' + subject + '_trial' + str(trial_num).zfill(3) + '_' + condition + '.csv'
			file = open(fileName, 'a')
			file.write(data_batch)
			file.close()		

	
			print 'End Trial ' + str(trial_num)
			
			# When trial_num is greater than TOTAL_TRIALS, end experiment
			if trial_num == TOTAL_TRIALS:
				is_freewalk = True
				goToStage('pretrial_00_01')
				trial_num = 1
			# Returns to Stage 1, resets clock
			else:
				trial_num += 1
				goToStage('pretrial_02_01')

	

# Restarts the loop, at a rate of 60Hz
viz.callback(viz.TIMER_EVENT,masterLoop)
viz.starttimer(0,1.0/90,viz.FOREVER)

