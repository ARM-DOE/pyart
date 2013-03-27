""" Cython wrapper around NASA TRMM RSL Library """

# data structures


cdef extern from "rsl.h":
    
    ctypedef unsigned short Range
    
    ctypedef struct Ray_header:
        int   month 
        int   day  
        int   year  
        int   hour  
        int   minute
        float sec 
        float unam_rng 
        float azimuth 
        int   ray_num   
        float elev    
        int   elev_num  
        int   range_bin1
        int   gate_size 
        float  vel_res  
        float sweep_rate
        int prf         
        int prf2        
        float azim_rate 
        float fix_angle 
        float pitch
        float roll      
        float heading   
        float pitch_rate
        float roll_rate 
        float heading_rate
        float lat
        float lon       
        int   alt      
        float rvc     
        float vel_east
        float vel_north
        float vel_up    
        int   pulse_count 
        float pulse_width 
        float beam_width 
        float frequency   
        float wavelength  
        float nyq_vel  
        float (*f)(Range x)     
        Range (*invf)(float x) 
        int   nbins

    ctypedef struct Ray:
        Ray_header h
        Range * range

    ctypedef struct Sweep_header:
        int sweep_num         
        float elev
        float azimuth 
        float beam_width 
        float vert_half_bw 
        float horz_half_bw
        int nrays
        float (*f)(Range x) 
        Range (*invf)(float x)

    ctypedef struct Sweep:           
        Sweep_header h   
        Ray **ray

    ctypedef struct Volume_header:
        char *type_str 
        int nsweeps
        float calibr_const   
        float (*f)(Range x)   
        Range (*invf)(float x) 

    ctypedef struct Volume:
        Volume_header h
        Sweep **sweep

    ctypedef struct Radar_header:
        int month, day, year
        int hour, minute
        float sec
        char radar_type[50]
        int nvolumes
        int number     
        char name[8]     
        char radar_name[8]
        char project[24] 
        char city[15]  
        char state[3]    
        char country[15]
        int latd
        int latm 
        int lats  
        int lond  
        int lonm  
        int lons   
        int height 
        int spulse 
        int lpulse
        int scan_mode 
        int vcp

    ctypedef struct Radar:
        Radar_header h
        Volume **v
    
    # Function definitions
    #Radar RSL_any_format_to_radar

    Radar * RSL_anyformat_to_radar(char *infile)
    void RSL_print_version()


