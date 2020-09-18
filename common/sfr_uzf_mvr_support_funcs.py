import numpy as np

# Defining a function that builds the new MF6 SFR connection information using
# original SFR input. This is a generalized function that can be used to
# convert MF2K5-based model to the new MF6 format.  Currently applied to the
# Sagehen and modsim models


def gen_mf6_sfr_connections(orig_seg, orig_rch):
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

            if rchx[4] == 1:  # checks if first rch of segment
                # Collect all segs that dump to the current one (there may not
                # be any)
                dumpersegs = [item for item in orig_seg if item[2] == segid]
                # For every seg that outflows to current, set last reach of it
                # as an upstream connection
                for dumper in dumpersegs:
                    dumper_seg_id = dumper[0]
                    rch_cnt = len(
                        [item for item in orig_rch if item[3] == dumper_seg_id]
                    )
                    lastrch = [
                        item
                        for item in orig_rch
                        if item[3] == dumper_seg_id and item[4] == rch_cnt
                    ]
                    idx = orig_rch.index(lastrch[0])
                    upconn.append(int(idx))

                # Current reach is the most upstream reach for current segment
                if iupseg == 0:
                    pass
                elif iupseg > 0:  # Lake connections, signified with negative
                    # numbers, aren't handled here
                    iupseg_rchs = [
                        item for item in orig_rch if item[3] == iupseg
                    ]
                    # Get the index of the last reach of the segement that was
                    # the upstream segment in the orig sfr file
                    idx = orig_rch.index(iupseg_rchs[len(iupseg_rchs) - 1])  #
                    upconn.append(idx)

                # Even if the first reach of a segement, it will have an outlet
                # either the next reach in the segment, or first reach of
                # outseg, which should be taken care of below
                if len(allrchs) > 1:
                    idx = orig_rch.index(rchx)
                    # adjust idx for 0-based and increment to next item in list
                    dnconn.append(int(idx + 1) * -1)

            elif rchx[4] > 1 and not rchx[4] == len(allrchs):
                # Current reach is 'interior' on the original segment and
                # therefore should only have 1 upstream & 1 downstream segement
                idx = orig_rch.index(rchx)
                # B/c 0-based, idx will already be incremented by -1
                upconn.append(int(idx - 1))
                # adjust idx for 0-based and increment to next item in list
                dnconn.append(int(idx + 1) * -1)  # all downstream connections
                # are negative in MF6

            if rchx[4] == len(allrchs):
                # If the last reach in a multi-reach segment, always need to
                # account for the reach immediately upstream (single reach segs
                # dealt with above), unless of course we're dealing with a
                # single reach segment like in the case of a spillway from a lk
                if len(allrchs) != 1:
                    idx = orig_rch.index(rchx)
                    # B/c 0-based, idx will already be incremented by -1
                    upconn.append(int(idx - 1))

                # Current reach is last reach in segment and may have multiple
                # downstream connections, particular when dealing with
                # diversions.
                if ioutseg == 0:
                    pass
                elif ioutseg > 0:  # Lake connections, signified with
                    # negative numbers, aren't handled here
                    idnseg_rchs = [
                        item
                        for item in orig_rch
                        if item[3] == ioutseg and item[4] == 1
                    ]
                    idx = orig_rch.index(idnseg_rchs[0])
                    # adjust idx for 0-based and increment to next item in list
                    dnconn.append(int(idx) * -1)

                # In addition to ioutseg, look for all segments that may have
                # the current segment as their iupseg
                possible_divs = [
                    item for item in orig_seg if item[3] == rchx[3]
                ]
                for segx in possible_divs:
                    # Next, peel out all first reach for any segments listed in
                    # possible_divs
                    first_rchs = [
                        item
                        for item in orig_rch
                        if item[3] == segx[0] and item[4] == 1
                    ]
                    for firstx in first_rchs:
                        idx = orig_rch.index(firstx)
                        # adjust idx for 0-based & increment to nxt itm in list
                        dnconn.append(int(idx) * -1)

            # Append the collection of upconn & dnconn as an entry in a list
            idx = orig_rch.index(rchx)
            # Adjust current index for 0-based
            conns.append([idx] + upconn + dnconn)

    return conns


def determine_runoff_conns_4mvr(pth, elev_arr, ibnd, orig_rch, nrow, ncol):

    # Get the sfr information stored in a companion script
    sfr_dat = orig_rch.copy()
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

            sfrlayout_conn_candidate_elev = 10000.0
            while True:
                direc = 0
                min_elev = elev_arr[curr_i, curr_j]

                # Look straight left
                if curr_j > 0:
                    if (
                        not sfrlayout[curr_i, curr_j - 1] == 0
                        and not ibnd[curr_i, curr_j - 1] == 0
                    ):  # Step in if neighbor is a stream cell
                        if elev_arr[curr_i, curr_j - 1] > 0 and (
                            elev_arr[curr_i, curr_j - 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i, curr_j - 1]
                            < sfrlayout_conn_candidate_elev
                        ):
                            sfrlayout_conn_candidate = sfrlayout[
                                curr_i, curr_j - 1
                            ]
                            sfrlayout_conn_candidate_elev = elev_arr[
                                curr_i, curr_j - 1
                            ]
                            stop_candidate = True

                    elif (
                        not elev_arr[curr_i, curr_j - 1] == 0
                        and not ibnd[curr_i, curr_j - 1] == 0
                    ):  # Step here if neighbor is not an sfr cell
                        if (
                            elev_arr[curr_i, curr_j - 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i, curr_j - 1] < min_elev
                        ):
                            elevcm1 = elev_arr[curr_i, curr_j - 1]
                            min_elev = elevcm1
                            direc = 2

                # Look up and left
                if curr_j > 0 and curr_i > 0:
                    if (
                        not sfrlayout[curr_i - 1, curr_j - 1] == 0
                        and not ibnd[curr_i - 1, curr_j - 1] == 0
                    ):  # Step in if neighbor is a stream cell
                        if elev_arr[curr_i - 1, curr_j - 1] > 0 and (
                            elev_arr[curr_i - 1, curr_j - 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i - 1, curr_j - 1]
                            < sfrlayout_conn_candidate_elev
                        ):
                            sfrlayout_conn_candidate = sfrlayout[
                                curr_i - 1, curr_j - 1
                            ]
                            sfrlayout_conn_candidate_elev = elev_arr[
                                curr_i - 1, curr_j - 1
                            ]
                            stop_candidate = True

                    elif (
                        not elev_arr[curr_i - 1, curr_j - 1] == 0
                        and not ibnd[curr_i - 1, curr_j - 1] == 0
                    ):  # Step here if neighbor is not an sfr cell
                        if (
                            elev_arr[curr_i - 1, curr_j - 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i - 1, curr_j - 1] < min_elev
                        ):
                            elevrm1cm1 = elev_arr[curr_i - 1, curr_j - 1]
                            min_elev = elevrm1cm1
                            direc = 5

                # Look straight right
                if curr_j < ncol - 1:
                    if (
                        not sfrlayout[curr_i, curr_j + 1] == 0
                        and not ibnd[curr_i, curr_j + 1] == 0
                    ):  # Step in if neighbor is a stream cell
                        if elev_arr[curr_i, curr_j + 1] > 0 and (
                            elev_arr[curr_i, curr_j + 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i, curr_j + 1]
                            < sfrlayout_conn_candidate_elev
                        ):
                            sfrlayout_conn_candidate = sfrlayout[
                                curr_i, curr_j + 1
                            ]
                            sfrlayout_conn_candidate_elev = elev_arr[
                                curr_i, curr_j + 1
                            ]
                            stop_candidate = True

                    elif (
                        not elev_arr[curr_i, curr_j + 1] == 0
                        and not ibnd[curr_i, curr_j + 1] == 0
                    ):  # Step here if neighbor is not an sfr cell
                        if (
                            elev_arr[curr_i, curr_j + 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i, curr_j + 1] < min_elev
                        ):
                            elevcm1 = elev_arr[curr_i, curr_j + 1]
                            min_elev = elevcm1
                            direc = 4

                # Look straight right and down
                if curr_i < nrow - 1 and curr_j < ncol - 1:
                    if (
                        not sfrlayout[curr_i + 1, curr_j + 1] == 0
                        and not ibnd[curr_i + 1, curr_j + 1] == 0
                    ):  # Step in if neighbor is a stream cell
                        if elev_arr[curr_i + 1, curr_j + 1] > 0 and (
                            elev_arr[curr_i + 1, curr_j + 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i + 1, curr_j + 1]
                            < sfrlayout_conn_candidate_elev
                        ):
                            sfrlayout_conn_candidate = sfrlayout[
                                curr_i + 1, curr_j + 1
                            ]
                            sfrlayout_conn_candidate_elev = elev_arr[
                                curr_i + 1, curr_j + 1
                            ]
                            stop_candidate = True

                    elif (
                        not elev_arr[curr_i + 1, curr_j + 1] == 0
                        and not ibnd[curr_i + 1, curr_j + 1] == 0
                    ):  # Step here if neighbor is not an sfr cell
                        if (
                            elev_arr[curr_i + 1, curr_j + 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i + 1, curr_j + 1] < min_elev
                        ):
                            elevrp1cp1 = elev_arr[curr_i + 1, curr_j + 1]
                            min_elev = elevrp1cp1
                            direc = 7

                # Look straight up
                if curr_i > 0:
                    if (
                        not sfrlayout[curr_i - 1, curr_j] == 0
                        and not ibnd[curr_i - 1, curr_j] == 0
                    ):  # Step in if neighbor is a stream cell
                        if elev_arr[curr_i - 1, curr_j] > 0 and (
                            elev_arr[curr_i - 1, curr_j]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i - 1, curr_j]
                            < sfrlayout_conn_candidate_elev
                        ):
                            sfrlayout_conn_candidate = sfrlayout[
                                curr_i - 1, curr_j
                            ]
                            sfrlayout_conn_candidate_elev = elev_arr[
                                curr_i - 1, curr_j
                            ]
                            stop_candidate = True

                    elif (
                        not elev_arr[curr_i - 1, curr_j] == 0
                        and not ibnd[curr_i - 1, curr_j] == 0
                    ):  # Step here if neighbor is not an sfr cell
                        if (
                            elev_arr[curr_i - 1, curr_j]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i - 1, curr_j] < min_elev
                        ):
                            elevcm1 = elev_arr[curr_i - 1, curr_j]
                            min_elev = elevcm1
                            direc = 3

                # Look up and right
                if curr_i > 0 and curr_j < ncol - 1:
                    if (
                        not sfrlayout[curr_i - 1, curr_j + 1] == 0
                        and not ibnd[curr_i - 1, curr_j + 1] == 0
                    ):  # Step in if neighbor is a stream cell
                        if elev_arr[curr_i - 1, curr_j + 1] > 0 and (
                            elev_arr[curr_i - 1, curr_j + 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i - 1, curr_j + 1]
                            < sfrlayout_conn_candidate_elev
                        ):
                            sfrlayout_conn_candidate = sfrlayout[
                                curr_i - 1, curr_j + 1
                            ]
                            sfrlayout_conn_candidate_elev = elev_arr[
                                curr_i - 1, curr_j + 1
                            ]
                            stop_candidate = True

                    elif (
                        not elev_arr[curr_i - 1, curr_j + 1] == 0
                        and not ibnd[curr_i - 1, curr_j + 1] == 0
                    ):  # Step here if neighbor is not an sfr cell
                        if (
                            elev_arr[curr_i - 1, curr_j + 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i - 1, curr_j + 1] < min_elev
                        ):
                            elevrm1cp1 = elev_arr[curr_i - 1, curr_j + 1]
                            min_elev = elevrm1cp1
                            direc = 6

                # Look straight down
                if curr_i < nrow - 1:
                    if (
                        not sfrlayout[curr_i + 1, curr_j] == 0
                        and not ibnd[curr_i + 1, curr_j] == 0
                    ):  # Step in if neighbor is a stream cell
                        if elev_arr[curr_i + 1, curr_j] > 0 and (
                            elev_arr[curr_i + 1, curr_j]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i + 1, curr_j]
                            < sfrlayout_conn_candidate_elev
                        ):
                            sfrlayout_conn_candidate = sfrlayout[
                                curr_i + 1, curr_j
                            ]
                            sfrlayout_conn_candidate_elev = elev_arr[
                                curr_i + 1, curr_j
                            ]
                            stop_candidate = True

                    elif (
                        not elev_arr[curr_i + 1, curr_j] == 0
                        and not ibnd[curr_i + 1, curr_j] == 0
                    ):  # Step here if neighbor is not an sfr cell
                        if (
                            elev_arr[curr_i + 1, curr_j]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i + 1, curr_j] < min_elev
                        ):
                            elevrp1 = elev_arr[curr_i + 1, curr_j]
                            min_elev = elevrp1
                            direc = 1

                # Look down and left
                if curr_i < nrow - 1 and curr_j > 0:
                    if (
                        not sfrlayout[curr_i + 1, curr_j - 1] == 0
                        and not ibnd[curr_i + 1, curr_j - 1] == 0
                    ):  # Step in if neighbor is a stream cell
                        if elev_arr[curr_i + 1, curr_j - 1] > 0 and (
                            elev_arr[curr_i + 1, curr_j - 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i + 1, curr_j - 1]
                            < sfrlayout_conn_candidate_elev
                        ):
                            sfrlayout_conn_candidate = sfrlayout[
                                curr_i + 1, curr_j - 1
                            ]
                            sfrlayout_conn_candidate_elev = elev_arr[
                                curr_i + 1, curr_j - 1
                            ]
                            stop_candidate = True

                    elif (
                        not elev_arr[curr_i + 1, curr_j - 1] == 0
                        and not ibnd[curr_i + 1, curr_j - 1] == 0
                    ):  # Step here if neighbor is not an sfr cell
                        if (
                            elev_arr[curr_i + 1, curr_j - 1]
                            < elev_arr[curr_i, curr_j]
                            and elev_arr[curr_i + 1, curr_j - 1] < min_elev
                        ):
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
                    # Check if encountered ibnd == 0, which may be a lake or
                    # boundary that drains out of model
                    if ibnd[curr_i, curr_j] == 0:
                        # This condition is dealt with after looping through
                        # all cells, see comment "Last step is set..."
                        break
                    pass  # Commence next downstream cell search

    # Last step is set the 0's in the vicinity of the lake equal to the
    # negative of the lake connection
    for i in np.arange(0, nrow):
        for j in np.arange(0, ncol):
            if sfrlayout_new[i, j] == 0 and ibnd[i, j] > 0:
                sfrlayout_new[i, j] = -1

    return sfrlayout_new
