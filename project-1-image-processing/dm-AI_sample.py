#**************************************************************************************

#

#   Driver Monitoring Systems using AI (code sample)

#

#   File: eyes_position.m

#   Author: Jacopo Sini

#   Company: Politecnico di Torino

#   Date: 19 Mar 2024

#

#**************************************************************************************



# 1 - Import the needed libraries

import cv2

import mediapipe as mp

import numpy as np 

import time

import statistics as st

import os

"COLORS"
blue =  (255, 0, 0)
green = (0, 255, 0)
red =   (0, 0, 255)
font_size = 0.8


# 2 - Set the desired setting

mp_face_mesh = mp.solutions.face_mesh

face_mesh = mp_face_mesh.FaceMesh(

    max_num_faces=1,

    refine_landmarks=True, # Enables  detailed eyes points

    min_detection_confidence=0.5,

    min_tracking_confidence=0.5
)

mp_drawing_styles = mp.solutions.drawing_styles

mp_drawing = mp.solutions.drawing_utils

drawing_spec = mp_drawing.DrawingSpec(thickness=1, circle_radius=1)

# Get the list of available capture devices (comment out)

if False:
    index = 0
    arr = []
    while True:
        dev = cv2.VideoCapture(index)
        try:
            arr.append(dev.getBackendName)
        except:
            break
        dev.release()
        index += 1
    print("* arr ", arr)


# 3 - Open the video source

cap = cv2.VideoCapture(0) # Local webcam (index start from 0)

# 4 - Iterate (within an infinite loop)

time_open_R = -1
time_open_L = -1
drowsy = -1

while cap.isOpened(): 

    # 4.1 - Get the new frame
    success, image = cap.read() 

    start = time.time()

    # Also convert the color space from BGR to RGB
    if image is None:
        break
        #continue
    #else: #needed with some cameras/video input format
        #image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

    # To improve performace
    image.flags.writeable = False

    # 4.2 - Run MediaPipe on the frame
    results = face_mesh.process(image)

    # To improve performance
    image.flags.writeable = True

    # Convert the color space from RGB to BGR
    #image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)

    img_h, img_w, img_c = image.shape

    point_RER = [] # Right Eye Right
    point_REB = [] # Right Eye Bottom
    point_REL = [] # Right Eye Left
    point_RET = [] # Right Eye Top

    point_LER = [] # Left Eye Right
    point_LEB = [] # Left Eye Bottom
    point_LEL = [] # Left Eye Left
    point_LET = [] # Left Eye Top

    point_REIC = [] # Right Eye Iris Center
    point_LEIC = [] # Left Eye Iris Center


    # 4.3 - Get the landmark coordinates

    if results.multi_face_landmarks:

        for face_landmarks in results.multi_face_landmarks:

            for idx, lm in enumerate(face_landmarks.landmark):

                # Eye Gaze (Iris Tracking)
                # Left eye indices list

                #LEFT_EYE =[ 362, 382, 381, 380, 374, 373, 390, 249, 263, 466, 388, 387, 386, 385,384, 398 ]

                # Right eye indices list

                #RIGHT_EYE=[ 33, 7, 163, 144, 145, 153, 154, 155, 133, 173, 157, 158, 159, 160, 161 , 246 ]

                #LEFT_IRIS = [473, 474, 475, 476, 477]

                #RIGHT_IRIS = [468, 469, 470, 471, 472]

                "R EYE CONTOUR"
                if idx == 33:
                    point_RER = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 145:
                    point_REB = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 133:
                    point_REL = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 159:
                    point_RET = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                "R EYE CONTOUR P2 P3 P5 P6"
                if idx == 160:
                    point_RP2 = (lm.x * img_w, lm.y * img_h)
                    cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 158:
                    point_RP3 = (lm.x * img_w, lm.y * img_h)
                    cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 153:
                    point_RP5 = (lm.x * img_w, lm.y * img_h)
                    cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 144:
                    point_RP6 = (lm.x * img_w, lm.y * img_h)
                    cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)



                "L EYE CONTOUR"
                if idx == 362:
                    point_LER = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 374:
                    point_LEB = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 263:
                    point_LEL = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 386:
                    point_LET = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)



                "L EYE CONTOUR P2 P3 P5 P6"
                if idx == 385:
                    point_LP2 = (lm.x * img_w, lm.y * img_h)
                    cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 387:
                    point_LP3 = (lm.x * img_w, lm.y * img_h)
                    cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 373:
                    point_LP5 = (lm.x * img_w, lm.y * img_h)
                    cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 380:
                    point_LP6 = (lm.x * img_w, lm.y * img_h)
                    cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)



                "R EYE PUPIL"
                if idx == 468:
                    point_REIC = (lm.x * img_w, lm.y * img_h)
                    cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=6, color=(255, 255, 0), thickness=-1)                    

                if idx == 469:
                    point_469 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=green, thickness=-1)

                if idx == 470:
                    point_470 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=green, thickness=-1)

                if idx == 471:
                    point_471 = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=green, thickness=-1)

                if idx == 472:
                    point_472 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=green, thickness=-1)



                "L EYE PUPIL"
                if idx == 473:
                    point_LEIC = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=(0, 255, 255), thickness=-1)

                if idx == 474:
                    point_474 = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=blue, thickness=-1)

                if idx == 475:
                    point_475 = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=blue, thickness=-1)

                if idx == 476:
                    point_476 = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=blue, thickness=-1)

                if idx == 477:
                    point_477 = (lm.x * img_w, lm.y * img_h)
                    #cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=blue, thickness=-1)





                # TODO
                if idx == 33 or idx == 263 or idx == 1 or idx == 61 or idx == 291 or idx == 199:
                    if idx == 1:
                        nose_2d = (lm.x * img_w, lm.y * img_h)
                        nose_3d = (lm.x * img_w, lm.y * img_h, lm.z * 3000)

                    x, y = int(lm.x * img_w), int(lm.y * img_h)

                # TODO

                #LEFT_IRIS = [473, 474, 475, 476, 477]
                if idx == 473 or idx == 362 or idx == 374 or idx == 263 or idx == 386: # iris points
                #if idx == 473 or idx == 474 or idx == 475 or idx == 476 or idx == 477: # eye border

                    if idx == 473:
                        left_pupil_2d = (lm.x * img_w, lm.y * img_h)
                        left_pupil_3d = (lm.x * img_w, lm.y * img_h, lm.z * 3000)

                    x, y = int(lm.x * img_w), int(lm.y * img_h)


                # TODO
                #RIGHT_IRIS = [468, 469, 470, 471, 472]
                if idx == 468 or idx == 33 or idx == 145 or idx == 133 or idx == 159: # iris points
                # if idx == 468 or idx == 469 or idx == 470 or idx == 471 or idx == 472: # eye border

                    if idx == 468:
                        right_pupil_2d = (lm.x * img_w, lm.y * img_h)
                        right_pupil_3d = (lm.x * img_w, lm.y * img_h, lm.z * 3000)

                    x, y = int(lm.x * img_w), int(lm.y * img_h)


            # 4.4. - Draw the positions on the frame

            l_eye_width = point_LEL[0] - point_LER[0]

            l_eye_height = point_LEB[1] - point_LET[1]

            l_eye_center = [(point_LEL[0] + point_LER[0])/2 ,(point_LEB[1] + point_LET[1])/2]

            horizontal_threshold = 0.1
            # cv2.circle(image, (int(l_eye_center[0]), int(l_eye_center[1])), radius=int(horizontal_threshold * l_eye_width), color=blue, thickness=-1) #center of eye and its radius 

            cv2.circle(image, (int(point_LEIC[0]), int(point_LEIC[1])), radius=3, color=green, thickness=-1) # Center of iris

            cv2.circle(image, (int(l_eye_center[0]), int(l_eye_center[1])), radius=2, color=(128, 128, 128), thickness=-1) # Center of eye

            # print("Left eye: x = " + str(np.round(point_LEIC[0],0)) + " , y = " + str(np.round(point_LEIC[1],0)))

            cv2.putText(image, "Left eye: x = " + str(np.round(point_LEIC[0],0)) + ", y = " + str(np.round(point_LEIC[1],0)), (200, 50), cv2.FONT_HERSHEY_SIMPLEX, font_size, green, 2) 



            r_eye_width = point_REL[0] - point_RER[0]

            r_eye_height = point_REB[1] - point_RET[1]

            r_eye_center = [(point_REL[0] + point_RER[0])/2 ,(point_REB[1] + point_RET[1])/2]

            #cv2.circle(image, (int(r_eye_center[0]), int(r_eye_center[1])), radius=int(horizontal_threshold * r_eye_width), color=blue, thickness=-1) #center of eye and its radius 

            cv2.circle(image, (int(point_REIC[0]), int(point_REIC[1])), radius=3, color=red, thickness=-1) # Center of iris

            cv2.circle(image, (int(r_eye_center[0]), int(r_eye_center[1])), radius=2, color=(128, 128, 128), thickness=-1) # Center of eye

            #print("right eye: x = " + str(np.round(point_REIC[0],0)) + " , y = " + str(np.round(point_REIC[1],0)))

            cv2.putText(image, "Right eye: x = " + str(np.round(point_REIC[0],0)) + ", y = " + str(np.round(point_REIC[1],0)), (200, 100), cv2.FONT_HERSHEY_SIMPLEX, font_size, red, 2) 

            # speed reduction (comment out for full speed)


            "EAR"
            R_EAR = (abs(point_RP2[1] - point_RP6[1]) + abs(point_RP3[1] - point_RP5[1])) / (2 * abs(point_REL[0] - point_RER[0]))
            L_EAR = (abs(point_LP2[1] - point_LP6[1]) + abs(point_LP3[1] - point_LP5[1])) / (2 * abs(point_LEL[0] - point_LER[0]))

            cv2.putText(image, "EAR: R = " + str(np.round(R_EAR, 3)) + ", L = " + str(np.round(L_EAR, 3)), (200, 150), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 


            "PERCLOS"
            EAR_max = 0.3

            if R_EAR > 0.8 * EAR_max:
                time_open_R += 1/15
            else:
                time_open_R = 0

            if L_EAR > 0.8 * EAR_max:
                time_open_L += 1/15
            else:
                time_open_L = 0

            if time_open_R > 2 or time_open_L > 2:
                drowsy = 1
            
            if drowsy > 0:
                cv2.putText(image, "DROWSY !!!", (250, 250), cv2.FONT_HERSHEY_SIMPLEX, 1, red, 2) 
                drowsy -= 1/15
            # cv2.putText(image, str(np.round(EAR_max, 3))+" "+str(np.round(time_open_R, 5))+" "+str(np.round(drowsy, 5)), (200, 200), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 
            # cv2.putText(image, str(np.round(EAR_max, 3))+" "+str(np.round(time_open_L, 5))+" "+str(np.round(drowsy, 5)), (200, 250), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 
        





            time.sleep(1/25) # [s]




        end = time.time()

        totalTime = end-start

        if totalTime>0:
            fps = 1 / totalTime
        else:
            fps=0

        #print("FPS:", fps)

        cv2.putText(image, f'FPS : {int(fps)}', (20,450), cv2.FONT_HERSHEY_SIMPLEX, 1.5, red, 2)

        # 4.5 - Show the frame to the user

        cv2.imshow('Technologies for Autonomous Vehicles - Driver Monitoring Systems using AI code sample', image)       

    if cv2.waitKey(5) & 0xFF == 27:
        break

# 5 - Close properly soruce and eventual log file

cap.release()

#log_file.close()

# [EOF]

