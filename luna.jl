# Code for the 2.4um, 65fs, 2.2mJ laser at UCF
#Two stages: first temporal compression, second tuning of an RDW emission in the UV/Vis

using Luna  #declare package to use propagation functions, need to add to Julia REPL before first use, see GitHub for documentation

# ##global variables:
    #parameters of input pulse using the Pharos laser
    λ0 = 2400e-9	    #central wavelength in meters
    τfwhm_initial= 65e-15	#pulse duration in seconds
    #parameters of the first fiber
    a1 = 268e-6		#core radius in meters
    flength1 = 1 	#fiber length in meters
    gas1 = :Kr		#gas in the fiber, Kr is krypton 
    energy1 = 500e-6	#pulse energy in Joules
    pressure1 = 1.0		#gas pressure in bar
    #parameters of the second fiber
    a2 = 175e-6		#core radius in meters
    flength2 = 1.0	#fiber length in meters
    compressor_stg1 =  Array{Any}(undef,1) #only using the first index but need array since mutuable in/out of functions
    compressor_stg2 =  Array{Any}(undef,6) #one index for each tested output pulse duration
    #plotting/grid parameters
    FTL=false #turns off the transformed limited case
    modes= 4 		 #total number of modes
    trange1 = 300e-15 	    #time grid size in seconds
    λrange1= (100e-9, 0.9λ0) #grid wavelength range in meters for plots for graphs
    λlims_stg1 = (100e-9, 1.1λ0) #grid wavelength range in meters for fibers for propagation
    trange2 = 500e-15  	    #time grid size in seconds 
    λlims_stg2 =(100e-9, 2λ0) #grid wavelength range in meters for fibers for propagation
    a=1 #arbituary array index to easily make variables mutuable for functions



#:'######:::'########:::::'###::::'########::'##::::'##::'######::
#'##... ##:: ##.... ##:::'## ##::: ##.... ##: ##:::: ##:'##... ##:
# ##:::..::: ##:::: ##::'##:. ##:: ##:::: ##: ##:::: ##: ##:::..::
# ##::'####: ########::'##:::. ##: ########:: #########:. ######::
# ##::: ##:: ##.. ##::: #########: ##.....::: ##.... ##::..... ##:
# ##::: ##:: ##::. ##:: ##.... ##: ##:::::::: ##:::: ##:'##::: ##:
#. ######::: ##:::. ##: ##:::: ##: ##:::::::: ##:::: ##:. ######::
#:......::::..:::::..::..:::::..::..:::::::::..:::::..:::......::

#graphs before compression 
function graphs_b4(func, range, bandpasss, TF) 
    Plotting.time_1D(func; modes=1,  FTL=FTL, trange=(-300e-15, 300e-15), bandpass=bandpasss) 		
        #^^plot of power as function of time for the mode specified 
         #can add bandpass= to look at pulse duration for certain wavelegnths
    Plotting.spec_1D(func; log10=TF, modes=1,  λrange=range)
        #^^ plots SED as function of wavelength for each mode
    Plotting.stats(func)
        #^^ plots lots of info about each mode
    Plotting.energy(func; modes=1, bandpass=bandpasss)
        #^^ graphs the energy as a function of distance for given wavelength range
end

#graphs after compression 
function graphs_after(func, propp, range, bandpasss, TF) #where propp corresponds to prop1, or prop2 for compression
    Plotting.time_1D(func; modes=1, trange=(-300e-15, 300e-15), FTL=FTL, propagate=propp, bandpass=bandpasss) 		
    #^^plot of power as function of time for the mode specified 
    #can add bandpass= to look at pulse duration for certain wavelegnths
    Plotting.spec_1D(func; log10=TF, modes=1, λrange=range)
    #^^ plots SED as function of wavelength for each mode
    Plotting.stats(func)
    #         #^^ plots lots of info about each mode
    Plotting.energy(func; modes=1, bandpass=bandpasss)
    #     #^^ graphs the energy as a function of distance for given wavelength range
end


## functions for optics, HCF1 and HCF2


#:::'##::::'######::'########:::::'######::'########::::'###:::::'######:::'########:
#:'####:::'##... ##:... ##..:::::'##... ##:... ##..::::'## ##:::'##... ##:: ##.....::
#:.. ##::: ##:::..::::: ##::::::: ##:::..::::: ##:::::'##:. ##:: ##:::..::: ##:::::::
#::: ##:::. ######::::: ##:::::::. ######::::: ##::::'##:::. ##: ##::'####: ######:::
#::: ##::::..... ##:::: ##::::::::..... ##:::: ##:::: #########: ##::: ##:: ##...::::
#::: ##:::'##::: ##:::: ##:::::::'##::: ##:::: ##:::: ##.... ##: ##::: ##:: ##:::::::
#:'######:. ######::::: ##:::::::. ######::::: ##:::: ##:::: ##:. ######::: ########:
#:......:::......::::::..:::::::::......::::::..:::::..:::::..:::......::::........::

#simulates the pulse down the first fiber
function stg1(gas, energy, pressure)
    pulse_in_stg1=Pulses.GaussPulse(;λ0, τfwhm=τfwhm_initial, energy=energy, power=nothing, ϕ=Float64[], m=1,mode=:lowest, polarisation=:linear, propagator=prop1!)
    #propagate through the first stage fiber to broaden bandwidth with the parameters listed above
    compressor_stg1[1]= prop_capillary(a1, flength1, gas, pressure; pulses=[pulse_in_stg1], λ0=λ0, trange=trange1, λlims=(100e-9, 1.1λ0), modes=modes)
end
        
function prop1!(Eω, grid) # this mutates its input as required for the LunaPulse
     Fields.prop_material!(Eω, grid, :SiO2, 1e-3, λ0) # 1 1-mm window
     Fields.prop_material!(Eω, grid, :Air, 1, λ0) # 2m of air
    Fields.prop_material!(Eω, grid, :KBr, 5e-3, λ0) # 1 1-mm window
   _, Eωopt = Fields.optcomp_material(Eω, grid, :SiO2, λ0, -2e-2, 2e-2) #wedges for optimizing compression, last numbers are  min_thickness then max_thickness
   Eω .= Eωopt
 end
    
# non-mutating version of prop! as required for plotting
function prop1(grid, Eω)
    Eωout = copy(Eω)
    prop1!(Eωout, grid)
    Eωout
end




#:'#######::'##::: ##:'########::::::'######::'########::::'###:::::'######:::'########:
#'##.... ##: ###:: ##: ##.... ##::::'##... ##:... ##..::::'## ##:::'##... ##:: ##.....::
#..::::: ##: ####: ##: ##:::: ##:::: ##:::..::::: ##:::::'##:. ##:: ##:::..::: ##:::::::
#:'#######:: ## ## ##: ##:::: ##::::. ######::::: ##::::'##:::. ##: ##::'####: ######:::
#'##:::::::: ##. ####: ##:::: ##:::::..... ##:::: ##:::: #########: ##::: ##:: ##...::::
# ##:::::::: ##:. ###: ##:::: ##::::'##::: ##:::: ##:::: ##.... ##: ##::: ##:: ##:::::::
# #########: ##::. ##: ########:::::. ######::::: ##:::: ##:::: ##:. ######::: ########:
#.........::..::::..::........:::::::......::::::..:::::..:::::..:::......::::........::

#simulates the pulse down the second fiber
function stg2(a, a2, gas, energy, pressure) 
    #obtain pulse from the first fiber end via LunaPulse
     pulse_in_stg2 = Pulses.LunaPulse(compressor_stg1[1]; energy, propagator=prop1!) #input pulse from the 1st stage
    #propagate through the second stage fiber to broaden bandwidth even more, with parameters listed above
    compressor_stg2[a]= prop_capillary(a2, flength2, gas, pressure; pulses=[pulse_in_stg2], λ0=λ0, trange=trange2, λlims=λlims_stg2, modes=modes)
end


#:::'##::::'######::'########:::::'######::'########::::'###:::::'######:::'########:
#:'####:::'##... ##:... ##..:::::'##... ##:... ##..::::'## ##:::'##... ##:: ##.....::
#:.. ##::: ##:::..::::: ##::::::: ##:::..::::: ##:::::'##:. ##:: ##:::..::: ##:::::::
#::: ##:::. ######::::: ##:::::::. ######::::: ##::::'##:::. ##: ##::'####: ######:::
#::: ##::::..... ##:::: ##::::::::..... ##:::: ##:::: #########: ##::: ##:: ##...::::
#::: ##:::'##::: ##:::: ##:::::::'##::: ##:::: ##:::: ##.... ##: ##::: ##:: ##:::::::
#:'######:. ######::::: ##:::::::. ######::::: ##:::: ##:::: ##:. ######::: ########:
#:......:::......::::::..:::::::::......::::::..:::::..:::::..:::......::::........::
bandpass=(0.5λ0, 1.5λ0) #for time graphs
λrange1=bandpass #for wavelength graphs

stg1(gas1, energy1, pressure1) #use predefined function to simulate pulse down first fiber


graphs_b4(compressor_stg1[1], λrange1, bandpass, false) #plot pulse after first fiber before compression- function, predefined λrange in global variables; last parameter sets log scale to true or false
Plotting.prop_2D(compressor_stg1[1], :λ, dBmin=-40.0,  λrange=λrange1, modes=1) #graphs color map of wavlength vs distance
graphs_after(compressor_stg1[1], prop1, λrange1, bandpass, false) #plot pulse after compression - prop1 for first fiber, predefined λrange in global variables 

#:'#######::'##::: ##:'########::::::'######::'########::::'###:::::'######:::'########:
#'##.... ##: ###:: ##: ##.... ##::::'##... ##:... ##..::::'## ##:::'##... ##:: ##.....::
#..::::: ##: ####: ##: ##:::: ##:::: ##:::..::::: ##:::::'##:. ##:: ##:::..::: ##:::::::
#:'#######:: ## ## ##: ##:::: ##::::. ######::::: ##::::'##:::. ##: ##::'####: ######:::
#'##:::::::: ##. ####: ##:::: ##:::::..... ##:::: ##:::: #########: ##::: ##:: ##...::::
# ##:::::::: ##:. ###: ##:::: ##::::'##::: ##:::: ##:::: ##.... ##: ##::: ##:: ##:::::::
# #########: ##::. ##: ########:::::. ######::::: ##:::: ##:::: ##:. ######::: ########:
#.........::..::::..::........:::::::......::::::..:::::..:::::..:::......::::........::

gas2=  :Ne #gas in the fiber, Ar is argon, He is Helium, Ne is Neon, Kr is Krypton
energy2= 350e-6 #in joules
pressure2= 9.0#in bar

stg2(a, a2, gas2, energy2, pressure2) #simulate second fiber


a=1

λrange2 = (150e-9, λ0) #for frequency graphs
bandpass2=λrange2 #for time graphs
graphs_b4(compressor_stg2[a], λrange2, bandpass2, false) #plot pulse after second fiber
Plotting.energy(compressor_stg2[a]; modes=1, bandpass=bandpass2)
Plotting.prop_2D(compressor_stg2[a], :λ, λrange=λrange2, modes=1) #graphs color map of wavlength vs distance

λrange2 = (240e-9, 370e-9) #for frequency graphs
bandpass2=λrange2 #for time graphs
graphs_b4(compressor_stg2[a], λrange2, bandpass2, false) #plot pulse after second fiber
Plotting.energy(compressor_stg2[a]; modes=1, bandpass=bandpass2)
Plotting.prop_2D(compressor_stg2[a], :λ, λrange=λrange2, modes=1) #graphs color map of wavlength vs distance


#scan range of energies for specfic pressure and output time/frequency graphs

pressure2=1.5 #set pressure value in bar
gas2=:Ar#gas
energy_scan1=120e-6 #beginning value of scan
energy_scan2=170e-6 #last value of scan
scan_steps=20e-6 #step to increment the scanned variables
x=energy_scan1 #set first value of for loop
a=2 #to differentiate the different energy scans
bandpass= (250e-9, 0.9λ0) #for time graphs
λrange1=bandpass #for wavelength graphs

for x in energy_scan1:scan_steps:energy_scan2
    stg2(a, a2, gas2, x, pressure2) #use predefined function to simulate pulse down first fiber
    graphs_b4(compressor_stg2[a], λrange1, bandpass, false) #output time and frequency graphs based on limits of λrange1 and bandpass
    Plotting.prop_2D(compressor_stg2[a], :λ, dBmin=-40.0,  λrange=λrange1, modes=1) #output color graph of wavelength vs distance
    a=a+1
end


