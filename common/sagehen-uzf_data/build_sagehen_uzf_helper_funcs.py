import os
import numpy as np
import pandas as pd

#             0, 1,  2,  3,    4,    5,    6,  7,     8,         9,     
orig_seg = [( 1, 1,  9,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg1"),
            ( 2, 1, 10,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg2"),
            ( 3, 1,  9,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg3"),
            ( 4, 1, 11,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg4"),
            ( 5, 1, 15,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg5"),
            ( 6, 1, 13,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg6"),
            ( 7, 1, 12,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg7"),
            ( 8, 1, 14,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg8"),
            ( 9, 1, 10,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg9"),
            (10, 1, 11,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg10"),
            (11, 1, 12,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg11"),
            (12, 1, 13,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg12"),
            (13, 1, 14,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg13"),
            (14, 1, 15,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg14"),
            (15, 1,  0,  0, 0.00,  0.0,  0.0,  0,  0.04, "origSeg15")]

orig_rch = [(1, 45,  9,  1,  1),
            (1, 44,  9,  1,  2),
            (1, 43,  9,  1,  3),
            (1, 43, 10,  1,  4),
            (1, 42, 10,  1,  5),
            (1, 42, 11,  1,  6),
            (1, 42, 12,  1,  7),
            (1, 41, 12,  1,  8),
            (1, 41, 13,  1,  9),
            (1, 40, 13,  1, 10),
            (1, 40, 14,  1, 11),
            (1, 39, 14,  1, 12),
            (1, 39, 15,  1, 13),
            (1, 39, 16,  1, 14),
            (1, 39, 17,  1, 15),
            (1, 38, 17,  1, 16),
            (1, 38, 18,  1, 17),
            (1, 38, 19,  1, 18),
            (1, 38, 20,  1, 19),
            (1, 38, 21,  1, 20),
            (1, 38, 22,  1, 21),
            (1, 38, 23,  1, 22),
            (1, 37, 24,  1, 23),
            (1, 37, 25,  1, 24),
            (1, 36, 25,  1, 25),
            (1, 36, 26,  1, 26),
            (1, 35, 26,  1, 27),
            (1, 35, 27,  1, 28),
            (1, 36, 28,  1, 29),
            (1, 36, 29,  1, 30),
            (1, 36, 30,  1, 31),
            (1, 36, 31,  1, 32),
            (1, 36, 32,  1, 33),
            (1, 35, 32,  1, 34),
            (1, 35, 33,  1, 35),
            (1, 34, 33,  1, 36),
            (1, 34, 34,  1, 37),
            (1, 33, 34,  1, 38),
            (1, 33, 35,  1, 39),
            (1, 32, 35,  1, 40),
            (1, 32, 36,  1, 41),
            (1, 32, 37,  1, 42),
            (1, 53, 37,  2,  1),
            (1, 52, 37,  2,  2),
            (1, 52, 38,  2,  3),
            (1, 51, 38,  2,  4),
            (1, 50, 38,  2,  5),
            (1, 50, 39,  2,  6),
            (1, 49, 39,  2,  7),
            (1, 48, 39,  2,  8),
            (1, 47, 39,  2,  9),
            (1, 47, 40,  2, 10),
            (1, 46, 40,  2, 11),
            (1, 45, 40,  2, 12),
            (1, 45, 39,  2, 13),
            (1, 44, 39,  2, 14),
            (1, 43, 38,  2, 15),
            (1, 42, 38,  2, 16),
            (1, 41, 38,  2, 17),
            (1, 40, 38,  2, 18),
            (1, 40, 39,  2, 19),
            (1, 39, 39,  2, 20),
            (1, 38, 39,  2, 21),
            (1, 37, 40,  2, 22),
            (1, 36, 40,  2, 23),
            (1, 36, 41,  2, 24),
            (1, 35, 41,  2, 25),
            (1, 34, 41,  2, 26),
            (1, 33, 41,  2, 27),
            (1, 32, 41,  2, 28),
            (1, 31, 42,  2, 29),
            (1, 31, 33,  3,  1),
            (1, 31, 34,  3,  2),
            (1, 31, 35,  3,  3),
            (1, 31, 36,  3,  4),
            (1, 31, 37,  3,  5),
            (1, 48, 48,  4,  1),
            (1, 47, 48,  4,  2),
            (1, 46, 48,  4,  3),
            (1, 46, 47,  4,  4),
            (1, 45, 47,  4,  5),
            (1, 44, 47,  4,  6),
            (1, 43, 47,  4,  7),
            (1, 42, 47,  4,  8),
            (1, 41, 47,  4,  9),
            (1, 41, 48,  4, 10),
            (1, 40, 48,  4, 11),
            (1, 39, 48,  4, 12),
            (1, 38, 47,  4, 13),
            (1, 37, 47,  4, 14),
            (1, 36, 48,  4, 15),
            (1, 35, 48,  4, 16),
            (1, 35, 49,  4, 17),
            (1, 34, 49,  4, 18),
            (1, 34, 50,  4, 19),
            (1, 33, 50,  4, 20),
            (1, 55, 72,  5,  1),
            (1, 54, 72,  5,  2),
            (1, 53, 72,  5,  3),
            (1, 52, 72,  5,  4),
            (1, 51, 72,  5,  5),
            (1, 50, 73,  5,  6),
            (1, 49, 73,  5,  7),
            (1, 48, 73,  5,  8),
            (1, 48, 74,  5,  9),
            (1, 47, 74,  5, 10),
            (1, 46, 75,  5, 11),
            (1, 45, 75,  5, 12),
            (1, 45, 76,  5, 13),
            (1, 44, 76,  5, 14),
            (1, 45, 62,  6,  1),
            (1, 44, 62,  6,  2),
            (1, 43, 62,  6,  3),
            (1, 43, 63,  6,  4),
            (1, 42, 63,  6,  5),
            (1, 41, 63,  6,  6),
            (1, 40, 63,  6,  7),
            (1, 24, 55,  7,  1),
            (1, 25, 55,  7,  2),
            (1, 25, 56,  7,  3),
            (1, 26, 56,  7,  4),
            (1, 27, 56,  7,  5),
            (1, 28, 57,  7,  6),
            (1, 29, 57,  7,  7),
            (1, 30, 57,  7,  8),
            (1, 31, 57,  7,  9),
            (1, 32, 57,  7, 10),
            (1, 33, 57,  7, 11),
            (1, 33, 58,  7, 12),
            (1, 34, 58,  7, 13),
            (1, 34, 59,  7, 14),
            (1, 35, 59,  7, 15),
            (1, 36, 59,  7, 16),
            (1, 37, 60,  7, 17),
            (1, 23, 71,  8,  1),
            (1, 24, 71,  8,  2),
            (1, 25, 71,  8,  3),
            (1, 26, 71,  8,  4),
            (1, 27, 72,  8,  5),
            (1, 27, 73,  8,  6),
            (1, 28, 73,  8,  7),
            (1, 29, 73,  8,  8),
            (1, 30, 73,  8,  9),
            (1, 31, 73,  8, 10),
            (1, 32, 73,  8, 11),
            (1, 33, 73,  8, 12),
            (1, 34, 73,  8, 13),
            (1, 34, 74,  8, 14),
            (1, 35, 74,  8, 15),
            (1, 36, 74,  8, 16),
            (1, 36, 73,  8, 17),
            (1, 37, 73,  8, 18),
            (1, 38, 72,  8, 19),
            (1, 39, 72,  8, 20),
            (1, 40, 72,  8, 21),
            (1, 41, 72,  8, 22),
            (1, 42, 72,  8, 23),
            (1, 42, 73,  8, 24),
            (1, 31, 38,  9,  1),
            (1, 31, 39,  9,  2),
            (1, 31, 40,  9,  3),
            (1, 31, 41,  9,  4),
            (1, 31, 42,  9,  5),
            (1, 30, 42,  9,  6),
            (1, 30, 43, 10,  1),
            (1, 30, 44, 10,  2),
            (1, 29, 44, 10,  3),
            (1, 29, 45, 10,  4),
            (1, 29, 46, 10,  5),
            (1, 29, 47, 10,  6),
            (1, 30, 47, 10,  7),
            (1, 30, 48, 10,  8),
            (1, 31, 49, 10,  9),
            (1, 32, 50, 10, 10),
            (1, 32, 51, 11,  1),
            (1, 33, 52, 11,  2),
            (1, 33, 53, 11,  3),
            (1, 34, 53, 11,  4),
            (1, 34, 54, 11,  5),
            (1, 35, 54, 11,  6),
            (1, 35, 55, 11,  7),
            (1, 35, 56, 11,  8),
            (1, 36, 57, 11,  9),
            (1, 36, 58, 11, 10),
            (1, 36, 59, 11, 11),
            (1, 37, 59, 11, 12),
            (1, 37, 60, 11, 13),
            (1, 38, 60, 11, 14),
            (1, 38, 61, 12,  1),
            (1, 38, 62, 12,  2),
            (1, 38, 63, 12,  3),
            (1, 39, 63, 12,  4),
            (1, 39, 64, 13,  1),
            (1, 39, 65, 13,  2),
            (1, 40, 65, 13,  3),
            (1, 40, 66, 13,  4),
            (1, 40, 67, 13,  5),
            (1, 40, 68, 13,  6),
            (1, 41, 69, 13,  7),
            (1, 41, 70, 13,  8),
            (1, 42, 71, 13,  9),
            (1, 42, 72, 13, 10),
            (1, 42, 73, 13, 11),
            (1, 42, 73, 14,  1),
            (1, 43, 73, 14,  2),
            (1, 43, 74, 14,  3),
            (1, 43, 75, 14,  4),
            (1, 44, 75, 14,  5),
            (1, 44, 76, 14,  6),
            (1, 44, 77, 15,  1),
            (1, 44, 78, 15,  2),
            (1, 44, 79, 15,  3),
            (1, 45, 79, 15,  4)
]

def get_sfrlist():
    return orig_rch

def gen_mf6_sfr_connections():
    conns = []
    for i in np.arange(0, len(orig_seg)):
        tup = orig_seg[i]
        segid = tup[0]
        ioutseg = tup[2]
        iupseg = tup[3]
        
        # Get all reaches associated with segment
        # Find an element in a list of tuples
        allrchs = [item for item in orig_rch if item[3] == segid]
        
        # Loop through allrchs and generate list of connections
        for rchx in allrchs:
            # rchx will be a tuple
            upconn = []
            dnconn = []
            
            if rchx[4] == 1:      # checks if first rch of segment
                # Collect all segs that dump to the current one (there may not be any)
                dumpersegs = [item for item in orig_seg if item[2] == segid]
                # For every seg that outflows to current, set last reach of it as
                # an upstream connection
                for dumper in dumpersegs:
                    dumper_seg_id = dumper[0]
                    rch_cnt = len([item for item in orig_rch if item[3] == dumper_seg_id])
                    lastrch = [item for item in orig_rch if item[3] == dumper_seg_id and item[4] == rch_cnt]
                    idx = orig_rch.index(lastrch[0])
                    upconn.append(int(idx))
                
                # Current reach is the most upstream reach for current segment
                if iupseg == 0:
                    pass
                elif iupseg > 0:  # Lake connections, signified with negative numbers, aren't handled here
                    iupseg_rchs = [item for item in orig_rch if item[3] == iupseg]
                    # Get the index of the last reach of the segement that was the upstream segment in the orig sfr file
                    idx = orig_rch.index(iupseg_rchs[len(iupseg_rchs)-1])   # From: https://stackoverflow.com/questions/20239312/find-an-exact-tuple-match-in-a-list-of-tuples-and-return-its-index
                    upconn.append(idx)
                
                # Even if the first reach of a segement, it will have an outlet,
                # either the next reach in the segment, or first reach of outseg, 
                # which should be taken care of below
                if len(allrchs) > 1:
                    idx = orig_rch.index(rchx)
                    # adjust idx for 0-based and increment to next item in list
                    dnconn.append(int(idx + 1) * -1)
            
            elif rchx[4] > 1 and not rchx[4] == len(allrchs):
                # Current reach is 'interior' on the original segment and therefore
                # should only have 1 upstream & 1 downstream segement
                idx = orig_rch.index(rchx)
                # B/c 0-based, idx will already be incremented by -1
                upconn.append(int(idx - 1))
                # adjust idx for 0-based and increment to next item in list
                dnconn.append(int(idx + 1) * -1)  # all downstream connections are negative in MF6
            
            if rchx[4] == len(allrchs):
                # If the last reach in a multi-reach segment, always need to account
                # for the reach immediately upstream (single reach segs dealt with 
                # above), unless of course we're dealing with a single reach segment
                # like in the case of a spillway from a lake
                if len(allrchs) != 1:
                    idx = orig_rch.index(rchx)
                    # B/c 0-based, idx will already be incremented by -1
                    upconn.append(int(idx - 1))
                
                # Current reach is last reach in segment and may have multiple 
                # downstream connections, particular when dealing with diversions.
                if ioutseg == 0:
                    pass
                elif ioutseg > 0:       # Lake connections, signified with negative numbers, aren't handled here
                    idnseg_rchs = [item for item in orig_rch if item[3] == ioutseg and item[4] == 1]
                    idx = orig_rch.index(idnseg_rchs[0])
                    # adjust idx for 0-based and increment to next item in list
                    dnconn.append(int(idx) * -1)
                    
                # In addition to ioutseg, look for all segments that may have the 
                # current segment as their iupseg
                possible_divs = [item for item in orig_seg if item[3] == rchx[3]]
                for segx in possible_divs:
                    # Next, peel out all first reach for any segments listed in possible_divs
                    first_rchs = [item for item in orig_rch if item[3] == segx[0] and item[4] == 1]
                    for firstx in first_rchs:
                        idx = orig_rch.index(firstx)
                        # adjust idx for 0-based and increment to next item in list
                        dnconn.append(int(idx) * -1)
            
            # Append the collection of upconn & dnconn as an entry in a list
            idx = orig_rch.index(rchx)
            # Adjust current index for 0-based
            conns.append([idx] + upconn + dnconn)
    
    return conns

def determine_runoff_conns_4mvr(pth, nrow, ncol):
    
    elev_arr = np.loadtxt(os.path.join(pth,'dis_support','top1.txt'), dtype=np.float)
    ibnd     = np.loadtxt(os.path.join(pth,'bas_support','ibnd1.txt'), dtype=np.int)
    
    # Get the sfr information stored in a companion script
    sfr_dat  = get_sfrlist()
    sfrlayout = np.zeros_like(ibnd)
    for i, rchx in enumerate(sfr_dat):
        row = rchx[1]
        col = rchx[2]
        sfrlayout[row - 1, col - 1] = i + 1
    
    sfrlayout_new = sfrlayout.copy()
    
    stop_candidate = False
    
    for i in np.arange(0, nrow):
        for j in np.arange(0, ncol):
        
            # Check to ensure current cell is active
            if ibnd[i, j] == 0:
                continue
            
            # Check to make sure it is not a stream cell
            if not sfrlayout[i, j] == 0:
                continue
            
            # Recursively trace path by steepest decent back to a stream
            curr_i = i
            curr_j = j
            
            sfrlayout_conn_candidate_elev = 10000.
            while True:
                direc = 0
                min_elev = elev_arr[curr_i, curr_j]
                
                # Look straight left
                if curr_j > 0:
                    if not sfrlayout[curr_i, curr_j - 1] == 0 and not ibnd[curr_i, curr_j - 1] == 0:   # Step in if neighbor is a stream cell
                        if elev_arr[curr_i, curr_j - 1] > 0 and (elev_arr[curr_i, curr_j - 1] < elev_arr[curr_i, curr_j] and
                                                                 elev_arr[curr_i, curr_j - 1] < sfrlayout_conn_candidate_elev):
                            sfrlayout_conn_candidate = sfrlayout[curr_i, curr_j - 1]
                            sfrlayout_conn_candidate_elev = elev_arr[curr_i, curr_j - 1]
                            stop_candidate = True
                    
                    elif not elev_arr[curr_i, curr_j - 1] == 0 and not ibnd[curr_i, curr_j - 1] == 0:  # Step here if neighbor is not an sfr cell
                        if elev_arr[curr_i, curr_j - 1] < elev_arr[curr_i, curr_j] and elev_arr[curr_i, curr_j - 1] < min_elev:
                            elevcm1 = elev_arr[curr_i, curr_j - 1]
                            min_elev = elevcm1
                            direc = 2
                
                # Look up and left
                if curr_j > 0 and curr_i > 0:
                    if not sfrlayout[curr_i - 1, curr_j - 1] == 0 and not ibnd[curr_i - 1, curr_j - 1] == 0:   # Step in if neighbor is a stream cell
                        if elev_arr[curr_i - 1, curr_j - 1] > 0 and (elev_arr[curr_i - 1, curr_j - 1] < elev_arr[curr_i, curr_j] and
                                                                     elev_arr[curr_i - 1, curr_j - 1] < sfrlayout_conn_candidate_elev):
                            sfrlayout_conn_candidate = sfrlayout[curr_i - 1, curr_j - 1]
                            sfrlayout_conn_candidate_elev = elev_arr[curr_i - 1, curr_j - 1]
                            stop_candidate = True

                    elif not elev_arr[curr_i - 1, curr_j - 1] == 0 and not ibnd[curr_i - 1, curr_j - 1] == 0:   # Step here if neighbor is not an sfr cell
                        if elev_arr[curr_i - 1, curr_j - 1] < elev_arr[curr_i, curr_j] and elev_arr[curr_i - 1, curr_j - 1] < min_elev:
                            elevrm1cm1 = elev_arr[curr_i - 1, curr_j - 1]
                            min_elev = elevrm1cm1
                            direc = 5

                
                # Look straight right
                if curr_j < ncol - 1:
                    if not sfrlayout[curr_i, curr_j + 1] == 0 and not ibnd[curr_i, curr_j + 1] == 0:   # Step in if neighbor is a stream cell
                        if elev_arr[curr_i, curr_j + 1] > 0 and (elev_arr[curr_i, curr_j + 1] < elev_arr[curr_i, curr_j] and
                                                                 elev_arr[curr_i, curr_j + 1] < sfrlayout_conn_candidate_elev):
                            sfrlayout_conn_candidate = sfrlayout[curr_i, curr_j + 1]
                            sfrlayout_conn_candidate_elev = elev_arr[curr_i, curr_j + 1]
                            stop_candidate = True
                    
                    elif not elev_arr[curr_i, curr_j + 1] == 0 and not ibnd[curr_i, curr_j + 1] == 0:  # Step here if neighbor is not an sfr cell
                        if elev_arr[curr_i, curr_j + 1] < elev_arr[curr_i, curr_j] and elev_arr[curr_i, curr_j + 1] < min_elev:
                            elevcm1 = elev_arr[curr_i, curr_j + 1]
                            min_elev = elevcm1
                            direc = 4
                
                # Look straight right and down
                if curr_i < nrow - 1 and curr_j < ncol - 1:
                    if not sfrlayout[curr_i + 1, curr_j + 1] == 0 and not ibnd[curr_i + 1, curr_j + 1] == 0:   # Step in if neighbor is a stream cell
                        if elev_arr[curr_i + 1, curr_j + 1] > 0 and (elev_arr[curr_i + 1, curr_j + 1] < elev_arr[curr_i, curr_j] and
                                                                     elev_arr[curr_i + 1, curr_j + 1] < sfrlayout_conn_candidate_elev):
                            sfrlayout_conn_candidate = sfrlayout[curr_i + 1, curr_j + 1]
                            sfrlayout_conn_candidate_elev = elev_arr[curr_i + 1, curr_j + 1]
                            stop_candidate = True
                    
                    elif not elev_arr[curr_i + 1, curr_j + 1] == 0 and not ibnd[curr_i + 1, curr_j + 1] == 0:   # Step here if neighbor is not an sfr cell
                        if elev_arr[curr_i + 1, curr_j + 1] < elev_arr[curr_i, curr_j] and elev_arr[curr_i + 1, curr_j + 1] < min_elev:
                            elevrp1cp1 = elev_arr[curr_i + 1, curr_j + 1]
                            min_elev = elevrp1cp1
                            direc = 7
                
                
                # Look straight up
                if curr_i > 0:
                    if not sfrlayout[curr_i - 1, curr_j] == 0 and not ibnd[curr_i - 1, curr_j] == 0:   # Step in if neighbor is a stream cell
                        if elev_arr[curr_i - 1, curr_j] > 0 and (elev_arr[curr_i - 1, curr_j] < elev_arr[curr_i, curr_j] and
                                                                 elev_arr[curr_i - 1, curr_j] < sfrlayout_conn_candidate_elev):
                            sfrlayout_conn_candidate = sfrlayout[curr_i - 1, curr_j]
                            sfrlayout_conn_candidate_elev = elev_arr[curr_i - 1, curr_j]
                            stop_candidate = True
                    
                    elif not elev_arr[curr_i - 1, curr_j] == 0 and not ibnd[curr_i - 1, curr_j] == 0:   # Step here if neighbor is not an sfr cell
                        if elev_arr[curr_i - 1, curr_j] < elev_arr[curr_i, curr_j] and elev_arr[curr_i - 1, curr_j] < min_elev:
                            elevcm1 = elev_arr[curr_i - 1, curr_j]
                            min_elev = elevcm1
                            direc = 3
                    
                
                # Look up and right
                if curr_i > 0 and curr_j < ncol - 1:
                    if not sfrlayout[curr_i - 1, curr_j + 1] == 0 and not ibnd[curr_i - 1, curr_j + 1] == 0:   # Step in if neighbor is a stream cell
                        if elev_arr[curr_i - 1, curr_j + 1] > 0 and (elev_arr[curr_i - 1, curr_j + 1] < elev_arr[curr_i, curr_j] and
                                                                 elev_arr[curr_i - 1, curr_j + 1] < sfrlayout_conn_candidate_elev):
                            sfrlayout_conn_candidate = sfrlayout[curr_i - 1, curr_j + 1]
                            sfrlayout_conn_candidate_elev = elev_arr[curr_i - 1, curr_j + 1]
                            stop_candidate = True
                    
                    elif not elev_arr[curr_i - 1, curr_j + 1] == 0 and not ibnd[curr_i - 1, curr_j + 1] == 0:   # Step here if neighbor is not an sfr cell
                        if elev_arr[curr_i - 1, curr_j + 1] < elev_arr[curr_i, curr_j] and elev_arr[curr_i - 1, curr_j + 1] < min_elev:
                            elevrm1cp1 = elev_arr[curr_i - 1, curr_j + 1]
                            min_elev = elevrm1cp1
                            direc = 6
                
                # Look straight down
                if curr_i < nrow - 1:
                    if not sfrlayout[curr_i + 1, curr_j] == 0 and not ibnd[curr_i + 1, curr_j] == 0:   # Step in if neighbor is a stream cell
                        if elev_arr[curr_i + 1, curr_j] > 0 and (elev_arr[curr_i + 1, curr_j] < elev_arr[curr_i, curr_j] and
                                                                 elev_arr[curr_i + 1, curr_j] < sfrlayout_conn_candidate_elev):
                            sfrlayout_conn_candidate = sfrlayout[curr_i + 1, curr_j]
                            sfrlayout_conn_candidate_elev = elev_arr[curr_i + 1, curr_j]
                            stop_candidate = True
                    
                    elif not elev_arr[curr_i + 1, curr_j] == 0 and not ibnd[curr_i + 1, curr_j] == 0:   # Step here if neighbor is not an sfr cell
                        if elev_arr[curr_i + 1, curr_j] < elev_arr[curr_i, curr_j] and elev_arr[curr_i + 1, curr_j] < min_elev:
                            elevrp1 = elev_arr[curr_i + 1, curr_j]
                            min_elev = elevrp1
                            direc = 1
                    
                # Look down and left
                if curr_i < nrow - 1 and curr_j > 0:
                    if not sfrlayout[curr_i + 1, curr_j - 1] == 0 and not ibnd[curr_i + 1, curr_j - 1] == 0:   # Step in if neighbor is a stream cell
                        if elev_arr[curr_i + 1, curr_j - 1] > 0 and (elev_arr[curr_i + 1, curr_j - 1] < elev_arr[curr_i, curr_j] and
                                                                 elev_arr[curr_i + 1, curr_j - 1] < sfrlayout_conn_candidate_elev):
                            sfrlayout_conn_candidate = sfrlayout[curr_i + 1, curr_j - 1]
                            sfrlayout_conn_candidate_elev = elev_arr[curr_i + 1, curr_j - 1]
                            stop_candidate = True
                    
                    elif not elev_arr[curr_i + 1, curr_j - 1] == 0 and not ibnd[curr_i + 1, curr_j - 1] == 0:   # Step here if neighbor is not an sfr cell
                        if elev_arr[curr_i + 1, curr_j - 1] < elev_arr[curr_i, curr_j] and elev_arr[curr_i + 1, curr_j - 1] < min_elev:
                            elevrp1cm1 = elev_arr[curr_i + 1, curr_j - 1]
                            min_elev = elevrp1cm1
                            direc = 8
                
                # if stop candidate found, don't move the cell indices
                if not stop_candidate:
                    # Direc corresponds to:
                    #  |----------------------
                    #  |  5  |    3    |  6  |
                    #  |----------------------
                    #  |  2  | cur_cel |  4  |
                    #  |----------------------
                    #  |  8  |    1    |  7  |
                    #  |----------------------
                    if direc == 0:
                        break
                    elif direc == 1:
                        curr_i += 1
                    elif direc == 2:
                        curr_j -= 1
                    elif direc == 3:
                        curr_i -= 1
                    elif direc == 4:
                        curr_j += 1
                    elif direc == 5:
                        curr_i -= 1
                        curr_j -= 1
                    elif direc == 6:
                        curr_i -= 1
                        curr_j += 1
                    elif direc == 7:
                        curr_i += 1
                        curr_j += 1
                    elif direc == 8:
                        curr_i += 1
                        curr_j -= 1
                
                if stop_candidate:
                    sfrlayout_new[i, j] = sfrlayout_conn_candidate
                    stop_candidate = False
                    break  # Bust out of while loop
                elif not stop_candidate:
                    # Check if encountered ibnd == 0, which may be a lake or boundary that drains out of model
                    if ibnd[curr_i, curr_j] == 0:
                        # This condition is dealt with after looping through all cells,
                        # see comment that starts, "Last step is set..."
                        break
                    pass  # Commence next downstream cell search

    # Last step is set the 0's in the vicinity of the lake equal to the negative of the lake connection
    for i in np.arange(0, nrow):
        for j in np.arange(0, ncol):
            if sfrlayout_new[i, j] == 0 and ibnd[i, j] > 0:
                sfrlayout_new[i, j] = -1
    
    # Once all cells are filled, save to array
    np.savetxt(os.path.join(pth, 'uzf_support','irunbnd_mf6.txt'), sfrlayout_new, fmt='%5d')
    
    return sfrlayout_new
