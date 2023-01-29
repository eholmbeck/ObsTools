import matplotlib.pyplot as plt
import numpy as np
import sys

from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time, TimeDelta
from timezonefinder import TimezoneFinder
from pytz import timezone

# ===============================================
# ===============================================

class Plan:
	def __init__(self, site='lco', year=2000, month=1, day=1, targets=[], exptimes=[]):
		self.site = self.get_site(site)
		self.year = year
		self.month= month
		self.day  = day
		
		self.start_local = None
		self.start_UT = None
		self.timezone_offset = None
		
		self.welcome(targets, exptimes)
		
	# ===========================================
	def welcome(self, targets, exptimes):
		sys.stdout.write('\033[2J\033[H')
		self.start()
		
		print('='*60)
		print(' '*15 + 'Welcome to my planner!')
		print('='*60 + '\n')
		print('One moment please...', end='\r')
		
		# Matolotlib is weird, so I'm not sure this is necessary to define here
		self.fig = plt.figure(figsize=(9,6))
		self.ax = self.fig.add_axes([0.08,0.1,0.55,0.8])
		plt.ion()
		self.initialize_plot()
		self.targets = []
		self.exptimes = []
		self.altitudes = []
		self.add_targets(targets)
		self.add_exptimes(exptimes)
				
		print('Target list and exposure times:')
		self.show_targets()
		
		print('\nOptions:')		
		self.show_options()
		
		self.start_time = self.dusk
		print('What would you like to do?')
		choice = input('')

		self.plan = []
		self.choose(choice)
		
	# ===========================================
	def refresh_screen(self, targets, exptimes):
		sys.stdout.write('\033[2J\033[H')
		
		print('='*60)
		print(' '*15 + 'Welcome to my planner!')
		print('='*60 + '\n')
		
		print('Target list and exposure times:')
		self.show_targets()
		
		print('\nOptions:')		
		self.show_options()
		
		print('What would you like to do?')
		choice = input('')
		self.choose(choice)
	
	# ===========================================
	def get_site(self, site):
		if type(site) is str:
			return EarthLocation.of_site(site)
		return site

	# ===========================================
	def find_timezone(self):
		tf = TimezoneFinder()
		return timezone(tf.timezone_at(lat=self.site.lat.value, lng=self.site.lon.value))

	# ===========================================
	def UT_to_local(self, ut_time):
		return ut_time - TimeDelta(self.timezone_offset.days*24*3600 + self.timezone_offset.seconds, format='sec')

	# ===========================================
	def local_to_UT(self, local_time):
		return local_time + TimeDelta(self.timezone_offset.days*24*3600 + self.timezone_offset.seconds, format='sec')
		
	# ===========================================
	# Sets some helpful variables
	def start(self, tres=100):
		self.tz = self.find_timezone()
	
		# Try at 6:00 pm in UT; find the offset, then reset the UT time
		run_start = {'year': self.year, 'month': self.month, 'day': self.day, 'hour': 18}
		time = Time(run_start, format='ymdhms', scale='utc')
		time_local = time.to_datetime(timezone=self.tz)
		self.timezone_offset = time_local.utcoffset()
		self.start_UT = self.UT_to_local(time) # Not a conversion, but a correction here
		self.start_local = self.start_UT.to_datetime(timezone=self.tz)
		
		self.dusk, self.dawn = self.twilight(-18)
		self.naut_dusk, self.naut_dawn = self.twilight(-12)
		
		self.times_jd = np.linspace(self.dusk.value, self.dawn.value, tres)
		night_length = (self.dawn - self.dusk).sec/3600.

		dusk_dt = self.dusk.to_datetime()
		t0 = dusk_dt.hour + dusk_dt.minute/60. + dusk_dt.second/3600.
		self.times_hr = np.linspace(t0, t0+night_length, tres)		
	
	# ===========================================
	def twilight(self, alt=-18, night=0):
		from scipy.optimize import fsolve
		x0 = (self.start_UT + TimeDelta(night*24*3600., format='sec')).jd
		dusk = fsolve(self.__find_twilight, x0=x0, args=alt)[0]
		dawn = fsolve(self.__find_twilight, x0=dusk+0.5, args=alt)[0]
		return Time([dusk, dawn], format='jd', scale='utc')

	# Time is in UT
	def __find_twilight(self, time, alt):
		time_object = Time(time, format='jd', scale='utc')
		altazframe = AltAz(obstime=time_object, location=self.site)
		sunaltaz = get_sun(time_object).transform_to(altazframe)
		return sunaltaz.alt.value - alt

	# ===========================================
	def sunset_sunrise(self, night=0):
		return self.twilight(alt=-2.25, night=night)
	
	# ===========================================
	# ===========================================
	def initialize_plot(self):
		offset = self.timezone_offset.days*24 + self.timezone_offset.seconds/3600.0
		
		self.ax.tick_params(labeltop=True, labelright=True)
		#thismanager = plt.get_current_fig_manager()
		#thismanager.window.wm_geometry("+900+0")
		
		xticks = np.arange(int(self.times_hr[0])-1, int(self.times_hr[-1])+2)
		self.ax.set_xlim(self.times_hr[0]-1, self.times_hr[-1]+1)
		self.ax.set_xticks(xticks)
		
		self.ax.set_ylim(0,90)
		self.ax.set_yticks(range(0,91,10))
		yticks = self.ax.get_yticks()
		'''
		# This would be so much cleaner if it worked..
		nxticks = len(xticks)
		for i in range(nxticks):
			self.ax.get_xticklabels()[i].set_x(xticks[i])
			self.ax.get_xticklabels()[i].set_text(r'${:.0f}$'.format(xticks[i]))
			self.ax.get_xticklabels()[i+nxticks].set_x(xticks[i])
			self.ax.get_xticklabels()[i+nxticks].set_y(1)
			self.ax.get_xticklabels()[i+nxticks].set_text(r'${:.0f}$'.format((xticks[i]+offset)%24))
		
		nyticks = len(yticks)
		for i in range(nyticks):
			self.ax.get_yticklabels()[i].set_y(yticks[i])
			self.ax.get_yticklabels()[i].set_text(r'${:.0f}^\circ$'.format(yticks[i]))
			self.ax.get_yticklabels()[i+nyticks].set_x(1)
			self.ax.get_yticklabels()[i+nyticks].set_y(yticks[i])
			self.ax.get_yticklabels()[i+nyticks].set_text([r'${:.2f}$'.format(1.0/np.sin(yticks[i]*np.pi/180.))])

		#self.ax.text(0.5, 1, r'Local Time ({:}, GMT${:+.0f})$'.format(self.tz.zone, offset),
		#		ha='center', va='top', transform=xlabel.get_transform())
		'''
		self.ax2 = self.ax.twinx()
		ax3 = self.ax.twiny()

		self.ax2.set_ylim(0,90)
		self.ax2.set_yticks(range(10, 91,10))
		self.ax2.set_yticklabels([r'${:.2f}$'.format(1.0/np.sin(y*np.pi/180.)) for y in self.ax2.get_yticks()])
		self.ax2.text(1.01,0,'Airmass', va='center', ha='left', rotation=90, transform=self.ax2.transAxes)

		ax3.set_xlabel(r'Local Time ({:}, GMT${:+.0f})$'.format(self.tz.zone,offset))
		ax3.set_xlim(self.ax.get_xlim())
		ax3.set_xticks(xticks)
		ax3.set_xticklabels([r'${:.0f}$'.format((ut_tick + offset)%24) for ut_tick in xticks])
		
		self.ax.grid(True, ls=':', color='gray')

		self.ax.set_ylabel('Altitude')
		self.ax2.text(1.01,0,'Airmass', va='center', ha='left', rotation=90, transform=self.ax2.transAxes)
		xlabel = self.ax.set_xlabel('UTC')
		
		self.ax2.axvline(self.times_hr[0], color='gray', alpha=0.6, ls='--')
		self.ax2.axvline(self.times_hr[-1], color='gray', alpha=0.6, ls='--')
		
		self.fig.canvas.draw()
		plt.pause(0.001)


	# ===========================================
	def add_track(self, target):
		coords = SkyCoord.from_name(target)
		altaz = [coords.transform_to(AltAz(obstime=Time(time, format='jd'),location=self.site))\
				 for time in self.times_jd]
		
		altitudes = np.array([alt.alt.value for alt in altaz])
		self.altitudes.append(altitudes)
		max_loc = np.where(altitudes==np.max(altitudes))[0][0]
		
		self.ax2.plot(self.times_hr, self.altitudes[-1], color='gray', lw=0.7, alpha=0.7)
		self.ax2.text(self.times_hr[max_loc], self.altitudes[-1][max_loc], len(self.targets), 
					  color='gray', fontsize='x-small', va='bottom')
		
		self.fig.canvas.draw()
		plt.pause(0.001)

	# ===========================================
	def add_targets(self, target_list):
		if type(target_list) is not list:
			target_list = list(target_list)
		for target in target_list:
			self.add_target(target)
	
	def add_target(self, single_target):
		self.targets.append(single_target)
		self.add_track(single_target)
	
	def add_exptimes(self, times_list):
		if type(times_list) is not list:
			times_list = list(times_list)
		for time in times_list:
			self.add_exptime(time)
	
	def add_exptime(self, single_time):
		self.exptimes.append(single_time)
	
	# ===========================================
	def show_targets(self):
		if len(self.targets)-len(self.exptimes) != 0:
			print('ERROR: Targets and exposure times are unequal lengths.')
			exit()
		
		longest_string = max([len(target) for target in self.targets])
		fmt = '{:>3.0f}: {:<%s} {:>6.0f} s ({:.2f} hr)' %(longest_string+2)
		for i in range(len(self.targets)):
			print(fmt.format(i+1, self.targets[i], self.exptimes[i], self.exptimes[i]/3600.))
		
		header = ' # | UTC   | Star_ID                  | Exp. (s)'
		fmt = '{:>2.0f} | {:02.0f}:{:02.0f} | {:<%s} |{:>8.0f}' %(longest_string+2)
		self.plan_text = [header]
		self.plan = []

	# ===========================================
	def show_options(self):
		fmt = '  {:<35}'
		print(fmt.format('[#]: Add target number # to plan and plot'))
		print(fmt.format('[p]: Print current plan')+fmt.format('[c]: Clear plan'))
		print(fmt.format('[u]: Undo')+fmt.format('[n]: Go to next night'))
		print(fmt.format('[q]: Quit')+fmt.format('[s]: Save plan to file\n'))

	# ===========================================
	def print_plan(self, nights):
		sys.stdout.write('\033[2A\033[1G\033[2K')
		sys.stdout.write('\033[%s;1H'%(13+len(self.targets)))
		sys.stdout.write('\033[0J')

		sunset, sunrise = self.sunset_sunrise().to_datetime()

		longest_string = max([len(target) for target in self.targets])		
		header_fmt = '{:<2} | {:5} | {:<%s} | {:<8}' %(longest_string+2)
		header = header_fmt.format('#', 'UT', 'Star ID', 'Exp. (s)')
		fmt = '{:>2.0f} | {:02.0f}:{:02.0f} | {:<%s} |{:>8.0f}' %(longest_string+2)
		
		print('Approximate Schedule:')
		print(header)
		print('-'*len(header))
		print(header_fmt.format('', '{:02.0f}:{:02.0f}'.format(sunset.hour, sunset.minute), 'Civil Sunset', ''))
		
		for p in self.plan:
			if len(p) != 0:
				print(fmt.format(*p))
			else:
				nights += 1
				print(header_fmt.format('', '{:02.0f}:{:02.0f}'.format(sunrise.hour, sunrise.minute), 'Civil Sunrise', ''))
				print('-'*len(header))
				sunset, sunrise = self.sunset_sunrise(night=nights).to_datetime()
				print(header_fmt.format('', '{:02.0f}:{:02.0f}'.format(sunset.hour, sunset.minute), 'Civil Sunset', ''))
		
		print(header_fmt.format('', '{:02.0f}:{:02.0f}'.format(sunrise.hour, sunrise.minute), 'Civil Sunrise', ''))
		print('-'*len(header))
		print('\nWhat would you like to do?\n')
		return nights
		
	# ===========================================
	# TODO: split this up
	def choose(self, choice):
		break_flag = False
		nights = 0
		while type(choice) is str:
			if choice.lower() not in 'ucns':
				try: choice = int(choice)-1
				except: break

				begin_exp = np.where(self.times_jd<=self.start_time.value)[0][-1]
				self.start_time += TimeDelta(self.exptimes[choice]+300.0, format='sec').jd # Include a 5-minute delay
				exptime = self.exptimes[choice]
				try: end_exp = np.where(self.times_jd>=self.start_time.value)[0][0]
				except:
					extra_time = (self.start_time-self.dawn).sec
					self.start_time -= TimeDelta(extra_time, format='sec')
					print('WARNING: {:} requires {:.0f}s more exposure past dawn!'.format(self.targets[choice], extra_time), end="\r")
					end_exp = len(self.times_jd)-1
					exptime -= extra_time
					break_flag = True
				
				self.ax.plot(self.times_hr[begin_exp:end_exp], self.altitudes[choice][begin_exp:end_exp],
							 label=self.targets[choice], color='C%i' %(len(self.plan)))
				self.ax.legend(loc='upper left', bbox_to_anchor=(1.1,1), handlelength=1.2, frameon=False)
	
				self.plan.append([choice+1, int(self.times_hr[begin_exp]), self.times_hr[begin_exp]%1 * 60,
							 	 self.targets[choice], exptime])
							
			elif choice in 'uU':
				try:
					self.start_time -= TimeDelta(self.plan[-1][-1]+300.0, format='sec')
					self.plan = self.plan[:-1]
				except IndexError:
					nights -= 1
					self.plan = self.plan[:-2]
					self.start_time -= TimeDelta(self.plan[-1][-1]+300.0, format='sec')

				# This needs to be done with some backends for some reason, otherwise the lines won't remove
				self.ax.lines[len(self.plan)-nights].set_data([],[])
				self.ax.lines.pop()

			elif choice in 'cC':
				for p in range(len(self.plan)-nights):
					self.ax.lines[-1].set_data([],[])
					self.ax.lines.pop()
				self.start_time = self.dusk
				self.plan = []
		
			elif choice in 'nN':
				self.plan.append([])
				self.start_time = self.dusk
				#nights = self.print_plan(nights)
			
			elif choice in 'sSaAeE':
				sys.stdout.write('\033[1A')
				print('This option has not yet been implemented :(\n')
			
			# Not yet implemented
			elif choice in 'sS':
				file_name = input('File name to save to: ')
			
			# Not working yet
			elif choice in 'aA':
				targets = input('Add targets either as a (Pythonic) list of strings or a single target string:\n')
				try: exec('targets=%s' %targets)
				except NameError: exec('targets=["%s"]' %targets)
				self.add_targets(targets)

				exptimes = input('Add exposure times (list or single):\n')
				try: exec('exptimes=%s' %exptimes)
				except NameError: exec('exptimes=[%s]' %exptimes)
				self.add_exptimes(exptimes)
					
				self.refresh_screen(targets, exptimes)
				
			# Not working yet
			elif choice in 'eE':
				targets = input('New target list (blank for same):\n')
				if targets != '':
					try: exec('targets=%s' %targets)
					except NameError: exec('targets=["%s"]' %targets)
					for p in self.targets:
						self.ax2.lines[-1].set_data([],[])
						self.ax2.lines.pop()
					self.targets = []
					self.add_targets(targets)
				
				self.exptimes = []
				exptimes = input('New exposure times:\n')
				try: exec('exptimes=%s' %exptimes)
				except NameError: exec('exptimes=[%s]' %exptimes)
				self.add_exptimes(exptimes)

				self.refresh_screen(targets, exptimes)

			else:
				break
	
			night = self.print_plan(nights)

			self.fig.canvas.draw()
			plt.pause(0.001)
	
			# Need to figure out a way to add this back in...
			#if break_flag:
			#	self.plan_text[-1] += '+'
			
			# see: https://gist.github.com/fnky/458719343aabd01cfb17a3a4f7296797
			sys.stdout.write('\033[1A\033[%sG\033[2K' %(' '*len(str(choice))))
			choice = input('')

