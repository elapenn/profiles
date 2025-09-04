#
#
import os
import sys
import argparse
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, timezone


# ---- SELECT SITE FROM THE LIST ----
SITES_list = ["Ispra EARLINET", "Potenza EARLINET"]
earID_list = ["ipr", "pot"]
monID_list = ["Ispra", "Potenza-EARLINET"]


# Set paths
plot_path = './graphs/'

#
def is_valid_date(date_string):
	# Chack string YYYY-MM-DD format
	date_obj = datetime.strptime(date_string, '%Y-%m-%d')

	if date_obj:
		return True
	else:
		return False


#
def is_valid_time(time_string):
    # Check string HH:MM format
	time_obj = datetime.strptime(time_string, '%H:%M')

	# Extract the hour
	hour = time_obj.hour

	# Check if the hour is in the allowed values
	allowed_hours = {0, 3, 6, 9, 12, 15, 18, 21} # UTC values
	if hour in allowed_hours:
		return True
	else:
		return False


#
def main(obs_site, obs_date, obs_hour):

	global plot_path

	global SITES_list
	global earID_list
	global monID_list

	#print(obs_site)

	# INDEX: 1 = Ispra; 2 = Potenza
	if obs_site == 'ISP':
		indx = 1
		print(SITES_list[indx-1])
	elif obs_site == 'POT':
		indx = 2
		print(SITES_list[indx-1])
	else:
		print('Only Potenza and Ispra sites allowed!')
		sys.exit(0)

	#
	USERS_SITE = SITES_list[indx-1]


	if not is_valid_date(obs_date):
		print('ERROR: Wrong date format!')
		sys.exit(0)

	if not is_valid_time(obs_hour):
		print('ERROR: Wrong hour format or values [UTC allowed: 0, 3, 6, 9, 12, 15, 18, 21]!')
		sys.exit(0)

	# Convert the string to a datetime object
	USERS_DT = datetime.strptime(str(obs_date+' '+obs_hour), '%Y-%m-%d %H:%M')
	#print('Chosen date: '+str(USERS_DT))
	#sys.exit(0)


	# MONARCH DATETIME
	#mon_file = "C:/Users/Admin/Desktop/ITINERIS/AERONET_ITALIA_25-06-2024/nc/MONARCH_Reanalysis_ec550du_" + monID_list[indx - 1] + "_2006-2016.nc"
	mon_file = "./DATA/MONARCH_Reanalysis_ec550du_" + monID_list[indx - 1] + "_2006-2016.nc"
	mon_DATETIME = nc.Dataset(mon_file).variables['DATETIME'][:]
	mon_DATETIME = np.array([datetime.fromtimestamp(dt, tz=timezone.utc).replace(tzinfo=None) for dt in mon_DATETIME])
	#print( mon_DATETIME )
	inmon_DT = np.where(mon_DATETIME == USERS_DT)[0]

	#print(mon_DATETIME, inmon_DT)
	#sys.exit()


	# EARLINET (IF SITE POT OR IPR)
	#ear_file = "C:/Users/Admin/Desktop/ITINERIS/Dust_Product_" + earID_list[indx - 1] + "_Lev02_DATETIME.nc"
	ear_file = "./DATA/Dust_Product_" + earID_list[indx - 1] + "_Lev02_DATETIME.nc"
	ear_DATETIMEsta = nc.Dataset(ear_file).variables['start_datetime'][:]
	#ear_DATETIMEsta = np.array([datetime.utcfromtimestamp(dt) for dt in ear_DATETIMEsta])
	ear_DATETIMEsta = np.array([datetime.fromtimestamp(dt, tz=timezone.utc).replace(tzinfo=None) for dt in ear_DATETIMEsta])

	ear_DATETIMEmid = nc.Dataset(ear_file).variables['middle_datetime'][:]
	#ear_DATETIMEmid = np.array([datetime.utcfromtimestamp(dt) for dt in ear_DATETIMEmid])
	ear_DATETIMEmid = np.array([datetime.fromtimestamp(dt, tz=timezone.utc).replace(tzinfo=None) for dt in ear_DATETIMEmid])

	ear_DATETIMEend = nc.Dataset(ear_file).variables['stop_datetime'][:]
	#ear_DATETIMEend = np.array([datetime.utcfromtimestamp(dt) for dt in ear_DATETIMEend])
	ear_DATETIMEend = np.array([datetime.fromtimestamp(dt, tz=timezone.utc).replace(tzinfo=None) for dt in ear_DATETIMEend])

	#print(ear_DATETIMEsta, ear_DATETIMEmid, ear_DATETIMEend)
	#sys.exit()

	dt_start = USERS_DT - timedelta(hours=3)
	dt_end = USERS_DT + timedelta(hours=3)
	inear_DT = np.where((ear_DATETIMEsta >= dt_start) & (ear_DATETIMEsta <= dt_end) |
						(ear_DATETIMEmid >= dt_start) & (ear_DATETIMEmid <= dt_end) |
						(ear_DATETIMEend >= dt_start) & (ear_DATETIMEend <= dt_end))[0]

	if len(inear_DT) > 1:
		closestIndex = np.argmin(np.abs(ear_DATETIMEmid[inear_DT] - USERS_DT))
		inear_DT = inear_DT[closestIndex]

	#inmon_DT = [] # Fake
	#inear_DT = [] # Fake


	# PLOT
	EB = 0  # 0: No errorbars; 1: With errorbars

	if len(inmon_DT) > 0 and len(inear_DT) == 0:  # CASE: MONARCH ONLY
		mon_ALTITUDE = nc.Dataset(mon_file).variables['ALTITUDE'][inmon_DT[0], :] # <-- MOD
		mon_DEC_AV = nc.Dataset(mon_file).variables['DEC_AV'][inmon_DT[0], :] # <-- MOD
		mon_DEC_STD = nc.Dataset(mon_file).variables['DEC_STD'][inmon_DT[0], :] # <-- MOD

		x_max = (np.max(mon_DEC_AV) + 0.07 * np.max(mon_DEC_AV)) * 1e3  # [km-1]

		plt.figure(1)
		plt.gcf().set_size_inches(6.76, 9.54)  # [width, height]

		if EB == 0:
			plt.plot(mon_DEC_AV * 1e3, mon_ALTITUDE / 1000, color=[.929, .694, .125], linewidth=4)
			plt.legend(['MONARCH'])
		else:
			plt.errorbar(mon_DEC_AV * 1e3, mon_ALTITUDE / 1000, xerr=mon_DEC_STD * 1e3, color=[.929, .694, .125], linewidth=4)
			plt.legend(['MONARCH'])

		plt.title(f"{USERS_SITE} {USERS_DT.strftime('%d-%b-%Y %H:%M')} UTC")
		plt.xlabel("[km$^{-1}$]")
		plt.ylabel('Height a.s.l. [km]')
		plt.grid(True)
		plt.ylim([0, 12])
		plt.xlim([0, x_max])
		plt.gca().tick_params(axis='both', which='major', labelsize=14)

		plt.savefig('my_plot.png', dpi=300)
		plt.show()
	#


	elif len(inmon_DT) == 0 and len(inear_DT) > 0:  # CASE: EARLINET ONLY
		l_EAR = np.loadtxt('./DATA/DATASET_L2_' + earID_list[indx - 1] + '532dust_start-202406_list.txt', dtype=str)
		ear_file = './DATA/desert_dust_prof_prod_collection_2015_2021/desert_dust_prof_prod_' + earID_list[indx - 1] + '_2015_2021/' + l_EAR[inear_DT][0]

		e_ALT = nc.Dataset(ear_file).variables['altitude'][:]
		e_DEC = nc.Dataset(ear_file).variables['dust_extinction'][0,0,:]
		e_DEC = np.nan_to_num(e_DEC)
		e_DECe = nc.Dataset(ear_file).variables['error_dust_extinction'][:]
		e_BOUNDS = nc.Dataset(ear_file).variables['time_bounds'][:]
		#e_BOUNDS = np.array([datetime.utcfromtimestamp(dt) for dt in e_BOUNDS])

		# Flatten the e_BOUNDS list and convert to a NumPy array
		e_B = np.array(e_BOUNDS).flatten()
		# Convert timestamps to datetime objects
		#e_BOUNDS = np.array([datetime.utcfromtimestamp(dt) for dt in e_B])
		e_BOUNDS = np.array([datetime.fromtimestamp(dt, tz=timezone.utc).replace(tzinfo=None) for dt in e_B])

		x_max = (np.max(e_DEC) + 0.07 * np.max(e_DEC)) * 1e3  # [km-1]

		plt.figure(2)
		plt.gcf().set_size_inches(6.76, 9.54)  # [width, height]

		if EB == 0:
			plt.plot(e_DEC * 1e3, e_ALT / 1000, color=[.85, .325, .098], linewidth=4)
			plt.legend([f'EARLINET {e_BOUNDS[0].strftime("%H:%M")}-{e_BOUNDS[1].strftime("%H:%M")}'])
		else:
			plt.errorbar(e_DEC * 1e3, e_ALT / 1000, xerr=e_DECe * 1e3, color=[.85, .325, .098], linewidth=4)
			plt.legend([f'EARLINET {e_BOUNDS[0].strftime("%H:%M")}-{e_BOUNDS[1].strftime("%H:%M")}'])

		plt.title(f"{USERS_SITE} {USERS_DT.strftime('%Y-%m-%d %H:%M')} UTC")
		plt.xlabel("[km$^{-1}$]")
		plt.ylabel('Height a.s.l. [km]')
		plt.grid(True)
		plt.ylim([0, 12])
		plt.xlim([0, x_max])
		plt.gca().tick_params(axis='both', which='major', labelsize=14)

		plt.savefig('my_plot.png', dpi=300)
		plt.show()
	#


	elif len(inmon_DT) > 0 and len(inear_DT) > 0:  # CASE: BOTH MONARCH & EARLINET

		# MONARCH
		mon_ALTITUDE = nc.Dataset(mon_file).variables['ALTITUDE'][inmon_DT[0], :]
		mon_DEC_AV = nc.Dataset(mon_file).variables['DEC_AV'][inmon_DT[0], :]
		mon_DEC_STD = nc.Dataset(mon_file).variables['DEC_STD'][inmon_DT[0], :]

		# EARLINET
		l_EAR = np.loadtxt('./DATA/DATASET_L2_' + earID_list[indx - 1] + '532dust_start-202406_list.txt', dtype=str)
		ear_file = './DATA/desert_dust_prof_prod_collection_2015_2021/desert_dust_prof_prod_' + earID_list[indx - 1] + '_2015_2021/' + l_EAR[inear_DT][0]

		e_ALT = nc.Dataset(ear_file).variables['altitude'][:]
		e_DEC = nc.Dataset(ear_file).variables['dust_extinction'][0,0,:]
		e_DEC = np.nan_to_num(e_DEC)
		e_DECe = nc.Dataset(ear_file).variables['error_dust_extinction'][:]
		e_BOUNDS = nc.Dataset(ear_file).variables['time_bounds'][:]

		# Flatten the e_BOUNDS list and convert to a NumPy array
		e_B = np.array(e_BOUNDS).flatten()
		# Convert timestamps to datetime objects
		#e_BOUNDS = np.array([datetime.utcfromtimestamp(dt) for dt in e_B])
		e_BOUNDS = np.array([datetime.fromtimestamp(dt, tz=timezone.utc).replace(tzinfo=None) for dt in e_B])

		x_max = (np.max([np.max(e_DEC), np.max(mon_DEC_AV)]) + 0.07 * np.max([np.max(e_DEC), np.max(mon_DEC_AV)])) * 1e3  # [km-1]

		plt.figure(3)
		plt.gcf().set_size_inches(6.76, 9.54)  # [width, height]

		if EB == 0:
			plt.plot(mon_DEC_AV * 1e3, mon_ALTITUDE / 1000, color=[.929, .694, .125], linewidth=4)
			plt.plot(e_DEC * 1e3, e_ALT / 1000, color=[.85, .325, .098], linewidth=4)
			plt.legend(['MONARCH', f'EARLINET {e_BOUNDS[0].strftime("%H:%M")}-{e_BOUNDS[1].strftime("%H:%M")}'])
		else:
			plt


		plt.title(f"{USERS_SITE} {USERS_DT.strftime('%Y-%m-%d %H:%M')} UTC")
		plt.xlabel("[km$^{-1}$]")
		plt.ylabel('Height a.s.l. [km]')
		plt.grid(True)
		plt.ylim([0, 12])
		plt.xlim([0, x_max])
		plt.gca().tick_params(axis='both', which='major', labelsize=14)

		output_file = plot_path+obs_site+'_'+obs_date.replace('-','')+'_'+obs_hour.replace(':','')+'.png'

		# Remove output file (if exists)
		if os.path.exists(output_file):
			os.remove(output_file)

		plt.savefig(output_file, dpi=300)
		#plt.show()
	#


#
def parse_args():

	"""
	Parses command-line arguments for...

	Returns:
		argparse.Namespace: Parsed command-line arguments
	"""
	parser = argparse.ArgumentParser(description="""
	Lidar profiles routine...
	Arguments:\n
		obs_site (str): Observation site to analyse.\n
		obs_date (str): Date in format YYYY-MM-DD.\n
		obs_hour (str): Hour in format HH:MM.\n

	Example usage:\n
		python3 main.py --obs_site POT --obs_date 2015-07-30 --obs_hour 21:00

	""")

	parser.add_argument('--obs_site', type=str, help='Observation site', required=True)
	parser.add_argument('--obs_date', type=str, help='Observation date', required=True)
	parser.add_argument('--obs_hour', type=str, help='Observation hour', required=True)

	parser.add_argument('--basepath', type=str, help='Base output path', default='.') # implicit

	if len(sys.argv) == 1:
		parser.print_help(sys.stderr)
		sys.exit(1)
	return parser.parse_args()


# Main program
if __name__ == '__main__':

	# Create the nested directories if they do not exist
	os.makedirs(plot_path, exist_ok=True)

	#
	args = parse_args()

	#
	print('###############################')
	print('#  Profiles plotting routine  #')
	print('###############################')
	#
	obs_site = args.obs_site
	obs_date = args.obs_date
	obs_hour = args.obs_hour

	#
	main(obs_site, obs_date, obs_hour)

#
