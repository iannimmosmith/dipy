""" Validating QuickBundles
"""

import numpy as np
import dipy as dp
# track reading
from dipy.io.dpy import Dpy
# segmenation
from dipy.segment.quickbundles import QuickBundles
# visualization
from fos import Window, Region
from fos.actor import Axes, Text3D
from bundle_picker import TrackLabeler
# metrics 
from dipy.tracking.metrics import downsample
from dipy.tracking.distances import (bundles_distances_mam,
					bundles_distances_mdf,
					most_similar_track_mam)

from matplotlib.mlab import find

def load_data(id):
	ids=['02','03','04','05','06','08','09','10','11','12']
	filename =  'data/subj_'+ids[id]+'_lsc_QA_ref.dpy'
	dp=Dpy(filename,'r')
	tracks=dp.read_tracks()
	dp.close()
	return tracks

def show_qb_streamlines(tracks,qb):
	# Create gui and message passing (events)
	w = Window(caption = 'QB validation', 
		width = 1200, 
		height = 800, 
		bgcolor = (0.,0.,0.2) )
	# Create a region of the world of actors
	region = Region( regionname = 'Main',
			extent_min = np.array([-5.0, -5, -5]),
			extent_max = np.array([5, 5 ,5]))
	# Create actors
	tl = TrackLabeler('Bundle Picker',
			qb,qb.downsampled_tracks(),
			vol_shape=(182,218,182),tracks_alpha=1)   
	ax = Axes(name = "3 axes", scale= 10, linewidth=2.0)
	vert = np.array( [[2.0,3.0,0.0]], dtype = np.float32 )
	ptr = np.array( [[.2,.2,.2]], dtype = np.float32 )
	tex = Text3D( "Text3D", vert, "(0,0,0)", 10*2.5, 10*.5, ptr)
	#Add actor to their region
	region.add_actor(ax)
	#region.add_actor(tex)
	region.add_actor(tl)
	#Add the region to the window
	w.add_region(region)
	w.refocus_camera()

	print 'Actors loaded'

	return w,region,ax,tex

def get_random_streamlines(tracks,N):	
	#qb = QuickBundles(tracks,dist,18)
	#N=qb.total_clusters()
	random_labels = np.random.permutation(np.arange(len(tracks)))[:N]
	random_streamlines = [tracks[i] for i in random_labels]
	return random_streamlines
		
def count_close_tracks(sla, slb, dist_thr=20):
        cnt_a_close = zeros(len(slb))
        for ta in sla:
            dta = bundles_distances_mdf([ta],slb)[0]
#            dta = bundles_distances_mam([ta],slb)[0]
            cnt_a_close += binarise(dta, dist_thr)
        return cnt_a_close

'''
coverage = # neighb tracks / #tracks 
         = cntT.sum()/len(T)

overlap = (cntT>1).sum()/len(T)

missed == (cntT==0).sum()/len(T)
'''

#virtuals/#tracks
        
'''
compare_streamline_sets(sla,slb,dist=20):
	d = bundles_distances_mdf(sla,slb)
	d[d<dist]=1
	d[d>=dist]=0
	return d 
'''

def binarise(D, thr):
#Replaces elements of D which are <thr with 1 and the rest with 0
        return 1*(np.array(D)<thr)

id=0

tracks=load_data(id)

track_subset_size = 50000

tracks=tracks[:track_subset_size]
print 'Streamlines loaded'
#qb=QuickBundles(tracks,20,18)
#print 'QuickBundles finished'
#print 'visualize/interact with streamlines'
#window,region,axes,labeler = show_qb_streamlines(tracks,qb)

qb = QuickBundles(tracks,20,18)
N=qb.total_clusters()
print 'QB finished with', N, 'clusters'

random_streamlines={}
for rep in [0]:
	random_streamlines[rep] = get_random_streamlines(qb.downsampled_tracks(), N)
	
# Thresholded distance matrices (subset x tracks) where subset Q = QB centroids
# and subset R = matched random subset. Matrices have 1 if the compared
# tracks have MDF distance < threshold a,d 0 otherwise.
#DQ=compare_streamline_sets(qb.virtuals(),qb.downsampled_tracks(), 20)
#DR=compare_streamline_sets(random_streamlines[0],qb.downsampled_tracks(), 20)

# The number of subset tracks 'close' to each track
#neighbours_Q = np.sum(DQ, axis=0)
#neighbours_R = np.sum(DR, axis=0)
neighbours_Q = count_close_tracks(qb.virtuals(), qb.downsampled_tracks(), 20)
neighbours_R = count_close_tracks(random_streamlines[0], qb.downsampled_tracks(), 20)

maxclose = np.int(np.max(np.hstack((neighbours_Q,neighbours_R))))

# The numbers of tracks 0, 1, 2, ... 'close' subset tracks
counts = [(np.int(n), len(find(neighbours_Q==n)), len(find(neighbours_R==n)))
          for n in range(maxclose+1)]

print np.array(counts)

# Typically counts_Q shows (a) very few tracks with 0 close QB
# centroids, (b) many tracks with a small number (between 1 and 3?) close QB
# tracks, and (c) few tracks with many (>3?) close QB tracks

# By contrast counts_R shows (a) a large number of tracks with 0 close
# R (random) neighbours, (b) fewer tracks with a small number of close R
# tracks, and (c) a long tail showing how the R sample has over-sampled
# in dense parts of the tractography, coming up with several rather
# similar tracks. By contast the QB tracks are dissimilar by design - or
# can be thought of as more evenly distributed in track space.

# The output below was generated with subject 02, 5k tracks, and threshold 20.
# Column 0 is the neighbour count, and Columns 1 and 2 are the
# number of tracks with that neighbour count.

# I suppose you could say this revealed some kind of sparseness for the
# QB subset by comparison with the Random one

#counts, Qfreq, Rfreq 
#191 clusters
#[[   0    2  668]
# [   1 1174  828]
# [   2 1905  720]
# [   3 1279  703]
# [   4  491  572]
# [   5  117  505]
# [   6   28  393]
# [   7    2  249]
# [   8    2  144]
# [   9    0  104]
# [  10    0   40]
# [  11    0   23]
# [  12    0   15]
# [  13    0   16]
# [  14    0    7]
# [  15    0    8]
# [  16    0    5]]

#The next table is for 10k tracks
#QB finished with 349 clusters
#[[   0    3  750]
# [   1 1777 1361]
# [   2 3674 1567]
# [   3 2928 1554]
# [   4 1195 1472]
# [   5  320 1130]
# [   6   80  800]
# [   7   19  643]
# [   8    3  375]
# [   9    1  192]
# [  10    0  112]
# [  11    0   35]
# [  12    0    7]
# [  13    0    2]]

#Here is a table for 50k tracks (took a long time!)
#QB finished with 480 clusters
#[[    0     7  2656]
# [    1  4726  4071]
# [    2 13934  4958]
# [    3 15927  5968]
# [    4  9912  6379]
# [    5  3901  5526]
# [    6  1204  4709]
# [    7   296  3798]
# [    8    55  3147]
# [    9    16  2815]
# [   10    21  2159]
# [   11     1  1479]
# [   12     0   966]
# [   13     0   621]
# [   14     0   337]
# [   15     0   208]
# [   16     0   121]
# [   17     0    61]
# [   18     0    16]
# [   19     0     3]
# [   20     0     2]]

# With 100k
#QB finished with 554 clusters
#[[    0    13  4260]
# [    1  7934  6270]
# [    2 23695  8319]
# [    3 31300 11154]
# [    4 22775 11716]
# [    5  9732 11578]
# [    6  3141  9930]
# [    7  1053  8880]
# [    8   262  7568]
# [    9    68  5511]
# [   10    26  3639]
# [   11     1  2838]
# [   12     0  2051]
# [   13     0  1936]
# [   14     0  1615]
# [   15     0  1102]
# [   16     0   713]
# [   17     0   461]
# [   18     0   242]
# [   19     0   119]
# [   20     0    60]
# [   21     0    32]
# [   22     0     6]]

'''
Now with thresholds (10,10) and 1k clusters
QB finished with 560 clusters
[[  0   0 197]
 [  1 874 388]
 [  2 112 229]
 [  3  13  94]
 [  4   1  55]
 [  5   0  15]
 [  6   0   8]
 [  7   0  14]]
'''

'''
Now with thresholds (10,10) and 10k tracks
QB finished with 1830 clusters
[[   0    6 1696]
 [   1 4879 1957]
 [   2 3454 1694]
 [   3 1315 1448]
 [   4  292 1060]
 [   5   52  687]
 [   6    2  401]
 [   7    0  285]
 [   8    0  200]
 [   9    0  164]
 [  10    0  114]
 [  11    0   80]
 [  12    0   43]
 [  13    0   39]
 [  14    0   27]
 [  15    0   21]
 [  16    0   20]
 [  17    0   19]
 [  18    0   21]
 [  19    0   12]
 [  20    0    8]
 [  21    0    3]
 [  22    0    1]]
'''

'''
and thresholds (10,10) and 100k tracks
QB finished with 3790 clusters
[[    0    39  8528]
 [    1 19543 10003]
 [    2 33795 10830]
 [    3 28215 10861]
 [    4 13257 10523]
 [    5  4121  9341]
 [    6   847  7964]
 [    7   156  6448]
 [    8    26  4819]
 [    9     1  3637]
 [   10     0  2997]
 [   11     0  2186]
 [   12     0  1891]
 [   13     0  1424]
 [   14     0  1291]
 [   15     0  1033]
 [   16     0   996]
 [   17     0   881]
 [   18     0   881]
 [   19     0   701]
 [   20     0   631]
 [   21     0   600]
 [   22     0   478]
 [   23     0   376]
 [   24     0   310]
 [   25     0   143]
 [   26     0   107]
 [   27     0    49]
 [   28     0    27]
 [   29     0    13]
 [   30     0     9]
 [   31     0     2]
 [   32     0     7]
 [   33     0     6]
 [   34     0     3]
 [   35     0     3]
 [   36     0     1]]
'''
