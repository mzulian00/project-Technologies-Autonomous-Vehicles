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
gray =  (128, 128, 128)
font_size = 0.6


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
    while True:
        dev = cv2.VideoCapture(index)
        try:
            arr.append(dev.getBackendName)
        except:
            break
        dev.release()
        index += 1


# 3 - Open the video source

cap = cv2.VideoCapture(0) # Local webcam (index start from 0)

# 4 - Iterate (within an infinite loop)

time_open_R = -1
time_open_L = -1
drowsy = -1

time_notfocus = 0
notfocus = -1

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

    face_2d = []
    face_3d = []

    right_eye_2d = []
    right_eye_3d = []

    left_eye_2d = []
    left_eye_3d = []

    # The camera matrix
    focal_length = 1 * img_w
    cam_matrix = np.array([ 
        [focal_length,  0,              img_h / 2],
        [0,             focal_length,   img_w / 2],
        [0,             0,              1        ]
    ])
    # The distorsion parameters
    dist_matrix = np.zeros((4, 1), dtype=np.float64)

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
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 158:
                    point_RP3 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 153:
                    point_RP5 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 144:
                    point_RP6 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)



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
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 387:
                    point_LP3 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 373:
                    point_LP5 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)

                if idx == 380:
                    point_LP6 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=red, thickness=-1)



                "RIGHT EYE PUPIL"
                if idx == 468:
                    point_REIC = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=(255, 255, 0), thickness=-1)                    

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



                "LEFT EYE PUPIL"
                if idx == 473:
                    point_LEIC = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=(0, 255, 255), thickness=-1)

                if idx == 474:
                    point_474 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=blue, thickness=-1)

                if idx == 475:
                    point_475 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=blue, thickness=-1)

                if idx == 476:
                    point_476 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=blue, thickness=-1)

                if idx == 477:
                    point_477 = (lm.x * img_w, lm.y * img_h)
                    # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=2, color=blue, thickness=-1)


                "FACE"
                if idx == 33 or idx == 263 or idx == 1 or idx == 61 or idx == 291 or idx == 199:
                    if idx == 1:
                        nose_2d = (lm.x * img_w, lm.y * img_h)
                        nose_3d = (lm.x * img_w, lm.y * img_h, lm.z * 3000)
                        # cv2.circle(image, (int(lm.x * img_w), int(lm.y * img_h)), radius=4, color=blue, thickness=-1)

                    x, y = int(lm.x * img_w), int(lm.y * img_h)

                    # Get the 2D Coordinates
                    face_2d.append([x, y])
                    # Get the 3D Coordinates
                    face_3d.append([x, y, lm.z])



                "LEFT IRIS"
                if idx == 473 or idx == 362 or idx == 374 or idx == 263 or idx == 386: # iris points
                    if idx == 473:
                        left_pupil_2d = (lm.x * img_w, lm.y * img_h)
                        left_pupil_3d = (lm.x * img_w, lm.y * img_h, lm.z * 3000)
                        cv2.circle(image, (int(left_pupil_2d[0]), int(left_pupil_2d[1])), radius=2, color=green, thickness=-1)

                    x, y = int(lm.x * img_w), int(lm.y * img_h)
                    left_eye_2d.append([x, y])
                    left_eye_3d.append([x, y, lm.z])


                "RIGHT IRIS"
                if idx == 468 or idx == 33 or idx == 145 or idx == 133 or idx == 159: # iris points
                    if idx == 468:
                        right_pupil_2d = (lm.x * img_w, lm.y * img_h)
                        right_pupil_3d = (lm.x * img_w, lm.y * img_h, lm.z * 3000)

                    x, y = int(lm.x * img_w), int(lm.y * img_h)
                    right_eye_2d.append([x, y])
                    right_eye_3d.append([x, y, lm.z])


            # 4.4. - Draw the positions on the frame

            l_eye_width = point_LEL[0] - point_LER[0]
            l_eye_height = point_LEB[1] - point_LET[1]
            l_eye_center = [(point_LEL[0] + point_LER[0])/2 ,(point_LEB[1] + point_LET[1])/2]

            horizontal_threshold = 0.1
            # cv2.circle(image, (int(l_eye_center[0]), int(l_eye_center[1])), radius=int(horizontal_threshold * l_eye_width), color=blue, thickness=-1) #center of eye and its radius 
            cv2.circle(image, (int(point_LEIC[0]), int(point_LEIC[1])), radius=3, color=green, thickness=-1) # Center of iris
            cv2.circle(image, (int(l_eye_center[0]), int(l_eye_center[1])), radius=2, color=gray, thickness=-1) # Center of eye
            # print("Left eye: x = " + str(np.round(point_LEIC[0],0)) + " , y = " + str(np.round(point_LEIC[1],0)))
            cv2.putText(image, "LEFT EYE:  x = " + str(np.round(point_LEIC[0],0)) + ", y = " + str(np.round(point_LEIC[1],0)), (270, 50), cv2.FONT_HERSHEY_SIMPLEX, font_size, green, 2) 



            r_eye_width = point_REL[0] - point_RER[0]
            r_eye_height = point_REB[1] - point_RET[1]
            r_eye_center = [(point_REL[0] + point_RER[0])/2 ,(point_REB[1] + point_RET[1])/2]

            #cv2.circle(image, (int(r_eye_center[0]), int(r_eye_center[1])), radius=int(horizontal_threshold * r_eye_width), color=blue, thickness=-1) #center of eye and its radius 
            cv2.circle(image, (int(point_REIC[0]), int(point_REIC[1])), radius=3, color=red, thickness=-1) # Center of iris
            cv2.circle(image, (int(r_eye_center[0]), int(r_eye_center[1])), radius=2, color=gray, thickness=-1) # Center of eye
            #print("right eye: x = " + str(np.round(point_REIC[0],0)) + " , y = " + str(np.round(point_REIC[1],0)))
            cv2.putText(image, "RIGHT EYE: x = " + str(np.round(point_REIC[0],0)) + ", y = " + str(np.round(point_REIC[1],0)), (270, 80), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2)



            "EAR"
            ear_R = (abs(point_RP2[1] - point_RP6[1]) + abs(point_RP3[1] - point_RP5[1])) / (2 * abs(point_REL[0] - point_RER[0]))
            ear_L = (abs(point_LP2[1] - point_LP6[1]) + abs(point_LP3[1] - point_LP5[1])) / (2 * abs(point_LEL[0] - point_LER[0]))
            # cv2.putText(image, "EAR: R = " + str(np.round(ear_R, 3)) + ", L = " + str(np.round(ear_L, 3)), (200, 150), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 



            "PERCLOS"
            EAR_max = 0.3

            if ear_R < 0.8 * EAR_max:
                time_open_R = time.time()

            if ear_L < 0.8 * EAR_max:
                time_open_L = time.time()

            time_threshold = 10
            if (time.time()-time_open_R > time_threshold or time.time()-time_open_L > time_threshold):
                drowsy = 1
            
            if drowsy > 0:
                cv2.putText(image, "DROWSY !", (20, 50), cv2.FONT_HERSHEY_SIMPLEX, font_size, red, 2) 
                drowsy -= 1/15
            else:
                cv2.putText(image, "AWAKE", (20, 50), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 

            # cv2.putText(image, str(np.round(EAR_max, 3))+" "+str(np.round(time_open_R, 5))+" "+str(np.round(drowsy, 5)), (200, 200), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 
            # cv2.putText(image, str(np.round(EAR_max, 3))+" "+str(np.round(time_open_L, 5))+" "+str(np.round(drowsy, 5)), (270, 120), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2)       


            "ROT MATRICES"
            face_2d = np.array(face_2d, dtype=np.float64)
            face_3d = np.array(face_3d, dtype=np.float64)
            left_eye_2d = np.array(left_eye_2d, dtype=np.float64)
            left_eye_3d = np.array(left_eye_3d, dtype=np.float64)
            right_eye_2d = np.array(right_eye_2d, dtype=np.float64)
            right_eye_3d = np.array(right_eye_3d, dtype=np.float64)

            # Solve PnP
            success, rot_vec, trans_vec = cv2.solvePnP(face_3d, face_2d, cam_matrix, dist_matrix)
            success_left_eye, rot_vec_left_eye, trans_vec_left_eye = cv2.solvePnP(left_eye_3d, left_eye_2d, cam_matrix, dist_matrix)
            success_right_eye, rot_vec_right_eye, trans_vec_right_eye = cv2.solvePnP(right_eye_3d, right_eye_2d, cam_matrix, dist_matrix)

            # Get rotational matrix
            rmat, jac = cv2.Rodrigues(rot_vec)
            rmat_left_eye, jac_left_eye = cv2.Rodrigues(rot_vec_left_eye)
            rmat_right_eye, jac_right_eye = cv2.Rodrigues(rot_vec_right_eye)

            # Get angles
            angles, mtxR, mtxQ, Qx, Qy, Qz = cv2.RQDecomp3x3(rmat)
            angles_left_eye, mtxR_left_eye, mtxQ_left_eye, Qx_left_eye, Qy_left_eye, Qz_left_eye = cv2.RQDecomp3x3(rmat_left_eye)
            angles_right_eye, mtxR_right_eye, mtxQ_right_eye, Qx_right_eye, Qy_right_eye, Qz_right_eye = cv2.RQDecomp3x3(rmat_right_eye)


            pitch = angles[0] * 1800
            yaw = -angles[1] * 1800
            roll = 180 + (np.arctan2(point_RER[1] - point_LEL[1], point_RER[0] - point_LEL[0]) * 180 / np.pi)
            if roll > 180:
                roll = roll - 360
            pitch_left_eye = angles_left_eye[0] * 1800
            yaw_left_eye = angles_left_eye[1] * 1800
            pitch_right_eye = angles_right_eye[0] * 1800
            yaw_right_eye = angles_right_eye[1] * 1800


            # Display directions
            lenght_line = 3

            nose_3d_projection, jacobian = cv2.projectPoints(nose_3d, rot_vec, trans_vec, cam_matrix, dist_matrix)
            p1 = (int(nose_2d[0]), int(nose_2d[1]))
            p2 = (int(nose_2d[0] - yaw * lenght_line), int(nose_2d[1] - pitch * lenght_line))
            cv2.line(image, p1, p2, gray, 2)

            left_eye_3d_projection, jacobian = cv2.projectPoints(left_pupil_3d, rot_vec_left_eye, trans_vec_left_eye, cam_matrix, dist_matrix)
            p1_left = (int(left_pupil_2d[0]), int(left_pupil_2d[1]))
            p2_left = (int(left_pupil_2d[0] - yaw * lenght_line), int(left_pupil_2d[1] - pitch * lenght_line))
            # cv2.line(image, p1_left, p2_left, blue, 3)

            right_eye_3d_projection, jacobian = cv2.projectPoints(right_pupil_3d, rot_vec_left_eye, trans_vec_left_eye, cam_matrix, dist_matrix)
            p1_right = (int(right_pupil_2d[0]), int(right_pupil_2d[1]))
            p2_right = (int(right_pupil_2d[0] - yaw * lenght_line), int(right_pupil_2d[1] - pitch * lenght_line))
            # cv2.line(image, p1_right, p2_right, blue, 3)


            "HEAD GAZING"
            # cv2.putText(image, "P = "+str(np.round(angles[0]*1800, 2))+" Y = "+str(np.round(-angles[1]*1800, 2))+" R = "+str(np.round(angles[2]*1800, 2)), (50, 300), cv2.FONT_HERSHEY_SIMPLEX, font_size/2, blue, 2) 
            p1 = (int(left_pupil_2d[0]), int(left_pupil_2d[1]))
            p2 = (int(right_pupil_2d[0]), int(right_pupil_2d[1]))
            # cv2.line(image, p1, p2, green, 1)


            if abs(angles[0]*1800) > 30 or abs(angles[1]*1800) > 30:
                cv2.putText(image, "HEAD NOT FOCUS!", (20, 80), cv2.FONT_HERSHEY_SIMPLEX, font_size, red, 2) 
            else:
                cv2.putText(image, "HEAD FOCUS", (20, 80), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 
            


            "EYE GAZING"
            cv2.line(image, p1_right, p2_right, gray, 2)
            diff_r = (point_REIC[0] - r_eye_center[0],  point_REIC[1] - r_eye_center[1])
            yaw, pitch = -diff_r[0]*5, diff_r[1]*5
            p2_right = (int(p2_right[0] - yaw * lenght_line), int(p2_right[1] - pitch * lenght_line))
            cv2.line(image, p1_right, p2_right, blue, 3)


            cv2.line(image, p1_left, p2_left, gray, 2)
            diff_l = (point_LEIC[0] - l_eye_center[0],  point_LEIC[1] - l_eye_center[1])
            yaw, pitch = -diff_l[0]*5, diff_l[1]*5
            p2_left = (int(p2_left[0] - yaw * lenght_line), int(p2_left[1] - pitch * lenght_line))
            cv2.line(image, p1_left, p2_left, green, 3)            

            val = 1.5
            if (np.linalg.norm(diff_r) < 3.5 or np.linalg.norm(diff_l) < 3.5) and (abs(diff_r[1]) < val and abs(diff_l[1]) < val):
                cv2.putText(image, "EYES FOCUS", (20, 110), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 
            else :
                cv2.putText(image, "EYES NOT FOCUS!", (20, 110), cv2.FONT_HERSHEY_SIMPLEX, font_size, red, 2) 
            
            # cv2.putText(image, str(np.round(diff_r[0], 3))+" "+str(np.round(diff_r[1], 5))+" "+ str(np.round(np.linalg.norm(diff_r), 5)), (250, 250), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 
            # cv2.putText(image, str(np.round(diff_l[0], 3))+" "+str(np.round(diff_l[1], 5))+" "+ str(np.round(np.linalg.norm(diff_l), 5)), (250, 290), cv2.FONT_HERSHEY_SIMPLEX, font_size, blue, 2) 




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

