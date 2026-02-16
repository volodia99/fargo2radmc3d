# import global variables
import par

import numpy as np
import os


# -------------------------
# script calling RADMC3D
# -------------------------
def write_radmc3d_script():
    
    # RT in the dust continuum
    if par.RTdust_or_gas == 'dust':
        command ='radmc3d image lambda '+str(par.wavelength*1e3)+' npix '+str(par.nbpixels)+' incl '+str(par.inclination)+' posang '+str(par.posangle+90.0)+' phi '+str(par.phiangle)
        if par.plot_tau == 'Yes':
            command ='radmc3d image tracetau lambda '+str(par.wavelength*1e3)+' npix '+str(par.nbpixels)+' incl '+str(par.inclination)+' posang '+str(par.posangle+90.0)+' phi '+str(par.phiangle)
        if par.polarized_scat == 'Yes':
            command=command+' stokes'

    # RT in gas lines
    if par.RTdust_or_gas == 'gas' or par.RTdust_or_gas == 'both':
        if par.widthkms == 0.0:
            command='radmc3d image iline '+str(par.iline)+' vkms '+str(par.vkms)+' npix '+str(par.nbpixels)+' incl '+str(par.inclination)+' posang '+str(par.posangle+90.0)+' phi '+str(par.phiangle)
        else:
            command='radmc3d image iline '+str(par.iline)+' widthkms '+str(par.widthkms)+' linenlam '+str(par.linenlam)+' npix '+str(par.nbpixels)+' incl '+str(par.inclination)+' posang '+str(par.posangle+90.0)+' phi '+str(par.phiangle)
            if par.plot_tau == 'Yes':
                command='radmc3d image tracetau iline '+str(par.iline)+' widthkms '+str(par.widthkms)+' linenlam '+str(par.linenlam)+' npix '+str(par.nbpixels)+' incl '+str(par.inclination)+' posang '+str(par.posangle+90.0)+' phi '+str(par.phiangle)
                #command='radmc3d tausurf 1.0 iline '+str(iline)+' widthkms '+str(widthkms)+' linenlam '+str(linenlam)+' npix '+str(nbpixels)+' incl '+str(inclination)+' posang '+str(posangle+90.0)+' phi '+str(phiangle)

    # optional: second-order ray tracing
    if par.secondorder == 'Yes':
        command=command+' secondorder'

    # write execution script
    if par.verbose == 'Yes':
        print(command)
    SCRIPT = open('script_radmc','w')
    '''
    if par.Tdust_eq_Thydro == 'No':
        SCRIPT.write('radmc3d mctherm; '+command)
    else:
        SCRIPT.write(command)        
    '''
    SCRIPT.write(command)
    SCRIPT.close()
    os.system('chmod a+x script_radmc')


# ---------------------------------------
# write spatial grid in file amr_grid.inp
# ---------------------------------------
def write_AMRgrid(F, R_Scaling=1, Plot=False):

    if par.verbose == 'Yes':
        print("writing spatial grid")
    path_grid='amr_grid.inp'

    grid=open(path_grid,'w')

    grid.write('1 \n')              # iformat/ format number = 1
    grid.write('0 \n')              # Grid style (regular = 0)
    grid.write('101 \n')            # coordsystem: 100 < spherical < 200 
    grid.write('0 \n')              # gridinfo
    grid.write('1 \t 1 \t 1 \n')    # incl x, incl y, incl z

    # spherical radius, colatitude, azimuth
    grid.write(str(F.nrad)+ '\t'+ str(F.ncol)+'\t'+ str(F.nsec)+'\n') 

    # nrad+1 dimension as we need to enter the coordinates of the cells edges
    for i in range(F.nrad + 1):  
        grid.write(str(F.redge[i]*F.culength*1e2)+'\t') # with unit conversion in cm
    grid.write('\n')

    # colatitude
    for i in range(F.ncol + 1):
        grid.write(str(F.tedge[i])+'\t')
    grid.write('\n')

    # azimuth
    for i in range(F.nsec + 1):
        grid.write(str(F.pedge[i])+'\t')
    grid.write('\n')
    
    grid.close()


# -----------------------
# writing out wavelength 
# -----------------------
def write_wavelength():
    wmin = 0.1
    wmax = 10000.0
    Nw = 150
    Pw = (wmax/wmin)**(1.0/(Nw-1))
    
    waves = np.zeros(Nw)
    waves[0] = wmin
    for i in range(1, Nw):
        waves[i]=wmin*Pw**i

    if par.verbose == 'Yes':
        print('writing wavelength_micron.inp')

    path = 'wavelength_micron.inp'
    wave = open(path,'w')
    wave.write(str(Nw)+'\n')
    for i in range(Nw):
        wave.write(str(waves[i])+'\n')
    wave.close()


# -----------------------
# writing out star parameters 
# -----------------------
def write_stars(Rstar = 1, Tstar = 6000):
    wmin = 0.1
    wmax = 10000.0
    Nw = 150
    Pw = (wmax/wmin)**(1.0/(Nw-1))
    
    waves = np.zeros(Nw)
    waves[0] = wmin
    for i in range(1, Nw):
        waves[i]=wmin*Pw**i

    if par.verbose == 'Yes':
        print('writing stars.inp')

    path = 'stars.inp'
    wave = open(path,'w')

    wave.write('\t 2\n')
    if par.central_binary == 'No':
        wave.write('1 \t'+str(Nw)+'\n')
        wave.write(str(Rstar*par.R_Sun)+'\t'+str(par.M_Sun)+'\t 0 \t 0 \t 0 \n')
    else:
        if par.fargo3d == 'Yes':
            import sys, subprocess

            command = par.awk_command+' " /^UNITOFLENGTHAU/ " '+par.dir+'/variables.par'
            # check which version of python we're using
            if sys.version_info[0] < 3:   # python 2.X
                buf = subprocess.check_output(command, shell=True)
            else:                         # python 3.X
                buf = subprocess.getoutput(command)
            culength_cm = float(buf.split()[1])*1.5e13  # from au to centimeters

            command = par.awk_command+' " /^UNITOFMASSMSUN/ " '+par.dir+'/variables.par'
            # check which version of python we're using
            if sys.version_info[0] < 3:   # python 2.X
                buf = subprocess.check_output(command, shell=True)
            else:                         # python 3.X
                buf = subprocess.getoutput(command)
            cumass_Msun = float(buf.split()[1])   # in solar masses

            # read planet0 and planet1.dat files which contain stars coordinates and mass
            f1, x_primary, y_primary, z_primary, f5, f6, f7, mass_primary, date, omega = np.loadtxt(par.dir+"/planet0.dat",unpack=True)
            f1, x_secondary, y_secondary, z_secondary, f5, f6, f7, mass_secondary, date, omega = np.loadtxt(par.dir+"/planet1.dat",unpack=True)

            # assume star radius proportionnel to M^1/3: R/Rsun = (M/Msun)^(1/3)
            r_primary   = (mass_primary/cumass_Msun)**(1./3)  # in solar radii
            r_secondary = (mass_secondary/cumass_Msun)**(1./3)  # in solar radii

            # finally write stars.inp with both stars
            wave.write('2 \t'+str(Nw)+'\n')
            # radius (cm), mass (g), x (cm), y (cm), z (cm)
            wave.write(str(r_primary[par.on]*par.R_Sun)+'\t'+str(mass_primary[par.on]*par.M_Sun/cumass_Msun)+'\t'+str(x_primary[par.on]*culength_cm)+'\t'+str(y_primary[par.on]*culength_cm)+'\t'+str(z_primary[par.on]*culength_cm)+'\n')
            wave.write(str(r_secondary[par.on]*par.R_Sun)+'\t'+str(mass_secondary[par.on]*par.M_Sun/cumass_Msun)+'\t'+str(x_secondary[par.on]*culength_cm)+'\t'+str(y_secondary[par.on]*culength_cm)+'\t'+str(z_secondary[par.on]*culength_cm)+'\n')


    for i in range(Nw):
        wave.write('\t'+str(waves[i])+'\n')
    wave.write('\t -'+str(Tstar)+'\n')
    if par.central_binary == 'Yes':
        wave.write('\t -'+str(Tstar)+'\n')
    wave.close()


# --------------------
# writing radmc3d.inp
# --------------------
def write_radmc3dinp(incl_dust = 1,
                     incl_lines = 0,
                     lines_mode = 1,
                     nphot = 1000000,
                     nphot_scat = 1000000,
                     nphot_spec = 1000000,
                     nphot_mono = 1000000,
                     istar_sphere = 0,
                     scattering_mode_max = 0,
                     tgas_eq_tdust = 1,
                     modified_random_walk = 0,
                     itempdecoup=1,
                     setthreads=2,
                     rto_style=3 ):

    if par.verbose == 'Yes':
        print('writing radmc3d.inp')

    RADMCINP = open('radmc3d.inp','w')
    inplines = ["incl_dust = "+str(int(incl_dust))+"\n",
                "incl_lines = "+str(int(incl_lines))+"\n",
                "lines_mode = "+str(int(lines_mode))+"\n",
                "nphot = "+str(int(nphot))+"\n",
                "nphot_scat = "+str(int(nphot_scat))+"\n",
                "nphot_spec = "+str(int(nphot_spec))+"\n",
                "nphot_mono = "+str(int(nphot_mono))+"\n",
                "istar_sphere = "+str(int(istar_sphere))+"\n",
                "scattering_mode_max = "+str(int(scattering_mode_max))+"\n",
                "tgas_eq_tdust = "+str(int(tgas_eq_tdust))+"\n",
                "modified_random_walk = "+str(int(modified_random_walk))+"\n",
                "itempdecoup = "+str(int(itempdecoup))+"\n",
                "setthreads="+str(int(setthreads))+"\n",
                "rto_style="+str(int(rto_style))+"\n"]

    RADMCINP.writelines(inplines)
    RADMCINP.close()

    
# --------------------
# writing lines.inp
# --------------------
def write_lines(specie,lines_mode):

    if par.verbose == 'Yes':
        print("writing lines.inp")
    path_lines='lines.inp'

    lines=open(path_lines,'w')

    lines.write('2 \n')              # <=== Put this to 2
    lines.write('1 \n')              # Nr of molecular or atomic species to be modeled
    # LTE calculations
    if lines_mode == 1:
        lines.write('%s    leiden    0    0    0'%specie)    # incl x, incl y, incl z
    else:
    # non-LTE calculations
        lines.write('%s    leiden    0    0    1\n'%specie)    # incl x, incl y, incl z
        lines.write('h2')
    lines.close()

    # Get molecular data file
    molecular_file = 'molecule_'+str(par.gasspecies)+'.inp'

    datafile = str(par.gasspecies)
    if par.gasspecies == 'hco+':
        datafile = 'hco+@xpol'
    if par.gasspecies == 'so':
        datafile = 'so@lique'
    if par.gasspecies == 'cs':
        datafile = 'cs'
    dat_file = datafile+'.dat'

    if (os.path.isfile(molecular_file) == False) and (os.path.isfile(dat_file) == False):

        # ---
        # check if curl is installed
        from shutil import which
        if which('curl') is None:
            sys.exit('curl is not installed on your system! I cannot download the molecular data file. Please install curl and restart!')
        # ---
    
        if par.verbose == 'Yes':
            print('--------- Downloading molecular data file ----------')
            
        command = 'curl -k -O https://home.strw.leidenuniv.nl/~moldata/datafiles/'+dat_file
        print(command)
        os.system(command)
        command = 'mv '+datafile+'.dat molecule_'+str(par.gasspecies)+'.inp'
        os.system(command)


# --------------------
# optional heating source due to viscous heating (heatsource.inp)
# --------------------
def write_heatsource_file():

    if par.hydro2D == 'No':
        gascube = par.gas.data*(par.gas.cumass*1e3)/((par.gas.culength*1e2)**3.)  # ncol, nrad, nsec, quantity is in g / cm^3
    else:
        gascube = par.gas.data*(par.gas.cumass*1e3)/((par.gas.culength*1e2)**2.)  # nrad, nsec, quantity is in g / cm^2
        # we now need to expand vertically, assuming a Gaussian distribution (copy paste of what is done is gas_density.py)

        # Allocate arrays
        rhogascube     = np.zeros((par.gas.ncol,par.gas.nrad,par.gas.nsec))
        rhogascube_cyl = np.zeros((par.gas.nver,par.gas.nrad,par.gas.nsec))

        # gas aspect ratio as function of r (or actually, R, cylindrical radius)
        hgas = par.aspectratio * (par.gas.rmed)**(par.flaringindex)
        hg2D = np.zeros((par.gas.nrad,par.gas.nsec))
        r2D  = np.zeros((par.gas.nrad,par.gas.nsec))
        for th in range(par.gas.nsec):
            hg2D[:,th] = hgas     # nrad, nsec
            r2D[:,th] = par.gas.rmed  # nrad, nsec

        # work out vertical expansion. First, for the array in cylindrical coordinates
        for j in range(par.gas.nver):
            rhogascube_cyl[j,:,:] = gascube * np.exp( -0.5*(par.gas.zmed[j]/hg2D/r2D)**2.0 )     # nver, nrad, nsec
            rhogascube_cyl[j,:,:] /= ( np.sqrt(2.*np.pi) * r2D * hg2D  * par.gas.culength*1e2)   # quantity is now in g cm^-3

        # then, sweep through the spherical grid 
        for j in range(par.gas.ncol):
            for i in range(par.gas.nrad):
                
                R = par.gas.rmed[i]*np.sin(par.gas.tmed[j])  # cylindrical radius
                icyl = np.argmin(np.abs(par.gas.rmed-R))
                if R < par.gas.rmed[icyl] and icyl > 0:
                    icyl-=1
                
                z = par.gas.rmed[i]*np.cos(par.gas.tmed[j])  # vertical altitude            
                jcyl = np.argmin(np.abs(par.gas.zmed-z))
                if z < par.gas.zmed[jcyl] and jcyl > 0:
                    jcyl-=1

                # bilinear interpolation
                if (icyl < par.gas.nrad-1 and jcyl < par.gas.nver-1 and icyl > 0):
                    dr = par.gas.rmed[icyl+1]-par.gas.rmed[icyl]
                    dz = par.gas.zmed[jcyl+1]-par.gas.zmed[jcyl]
                    
                    xij     = (par.gas.rmed[icyl+1]-R) * (par.gas.zmed[jcyl+1]-z) / (dr*dz)
                    if xij < 0 or xij > 1:
                        print('beware that xij < 0 or xij > 1:',i,j,xij,par.gas.rmed[icyl+1]-R,dr,par.gas.zmed[jcyl+1]-z,dz)
                    
                    xijp1   = (par.gas.rmed[icyl+1]-R) * (z-par.gas.zmed[jcyl])   / (dr*dz)
                    if xijp1 < 0 or xijp1 > 1:
                        print('beware that xijp1 < 0 or xijp1 > 1:',i,j,xijp1,par.gas.rmed[icyl+1]-R,dr,z-par.gas.zmed[jcyl],dz)

                    xip1j   = (R-par.gas.rmed[icyl])   * (par.gas.zmed[jcyl+1]-z) / (dr*dz)
                    if xip1j < 0 or xip1j > 1:
                        print('beware that xip1j < 0 or xip1j > 1:',i,j,xip1j,R-par.gas.rmed[icyl],dr,par.gas.zmed[jcyl+1]-z,dz)

                    xip1jp1 = (R-par.gas.rmed[icyl])   * (z-par.gas.zmed[jcyl])   / (dr*dz)
                    if xip1jp1 < 0 or xip1jp1 > 1:
                        print('beware that xip1jp1 < 0 or xip1jp1 > 1:',i,j,xip1jp1,R-par.gas.rmed[icyl],dr,z-par.gas.zmed[jcyl],dz)

                    rhogascube[j,i,:] = rhogascube_cyl[jcyl,icyl,:]*xij +\
                    rhogascube_cyl[jcyl+1,icyl,:]*xijp1 +\
                    rhogascube_cyl[jcyl,icyl+1,:]*xip1j +\
                    rhogascube_cyl[jcyl+1,icyl+1,:]*xip1jp1
                
                else:
                    # simple nearest-grid point interpolation...
                    rhogascube[j,i,:] = rhogascube_cyl[jcyl,icyl,:]

    HSRC = open('heatsource.binp','wb')
   # requested header
    # hdr[0] = format number
    # hdr[1] = data precision (8 means double)
    # hdr[2] = nb of grid cells
    hdr = np.array([1, 8, par.gas.nrad*par.gas.nsec*par.gas.ncol], dtype=int)
    hdr.tofile(HSRC)

    # Default case: uniform microturbulence set by 'turbvel' parameter in params.dat
    visc_heating_rate = np.zeros((par.gas.ncol,par.gas.nrad,par.gas.nsec))  # ncol, nrad, nsec in erg/cm^3/s

    # model used in circumbinary discs simulations (2026)
    for i in range(par.gas.nrad):
        if par.gas.rmed[i] < 3.5:
            myalpha = 0.05 # 0.05   # inside cavity
        else:
            myalpha = 0.001 # 1e-3   # outside cavity
        for j in range(par.gas.ncol):
            # cylindrical radius
            r = par.gas.rmed[i] * np.sin(par.gas.tmed[j])
            # aspect ratio
            h = par.aspectratio * r**(par.flaringindex)
            # Keplerian velocity in cm/s
            vk = np.sqrt(par.G * (par.gas.cumass*1e3) / (r*1e2*par.gas.culength))
            # viscous heating rate = nu rho (r d_r Omega)^2 = 9/4 nu rho Omega^2 with nu = alpha cs^2 / Omega
            visc_heating_rate[j,i,:] = (9./4) * myalpha * rhogascube[j,i,:] * (h**2.0) * (vk**3.0) / (r*1e2*par.gas.culength)
            # if j == par.gas.ncol//2-1:
            #     print(visc_heating_rate[j,i,0],r,myalpha,h,vk,rhogascube[j,i,0])

    # If writing data in an ascii file the ordering should be: nsec, ncol, nrad.
    # We therefore need to swap axes of array visc_heating_rate
    # before dumping it in a binary file! just like mastermind game!
    visc_heating_rate = np.swapaxes(visc_heating_rate, 0, 1)  # nrad ncol nsec
    visc_heating_rate = np.swapaxes(visc_heating_rate, 0, 2)  # nsec ncol nrad
    visc_heating_rate.tofile(HSRC)
    HSRC.close()