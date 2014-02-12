""" Cython Declarations for NASA TRMM RSL rsl.h file """

cdef extern from "rsl.h":
   
    cdef char[6] _RSL_VERSION_STR "RSL_VERSION_STR"

    # data structures
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
    
    # function definitions
    Radar * RSL_anyformat_to_radar(char *infile)
    Radar * RSL_wsr88d_to_radar(char *infile, char *call_or_first_tape_file)
    
    Volume * RSL_copy_volume(Volume *v)
    
    Volume * RSL_new_volume(int max_sweeps)
    Sweep * RSL_new_sweep(int max_rays) 
    Ray * RSL_new_ray(int max_bins)
   
    void RSL_free_ray(Ray *r)
    void RSL_free_sweep(Sweep *s)
    void RSL_free_volume(Volume *v)
    void RSL_free_radar(Radar *r)
    
    void RSL_print_version()
    
    # conversion functions
    float DZ_F(Range x)
    float VR_F(Range x)
    float SW_F(Range x)
    float CZ_F(Range x)
    float ZT_F(Range x)
    float DR_F(Range x)
    float LR_F(Range x)
    float ZD_F(Range x)
    float DM_F(Range x)
    float RH_F(Range x)
    float PH_F(Range x)
    float XZ_F(Range x)
    float CD_F(Range x)
    float MZ_F(Range x)
    float MD_F(Range x)
    float ZE_F(Range x)
    float VE_F(Range x)
    float KD_F(Range x)
    float TI_F(Range x)
    float DX_F(Range x)
    float CH_F(Range x)
    float AH_F(Range x)
    float CV_F(Range x)
    float AV_F(Range x)
    float SQ_F(Range x)
    float VS_F(Range x)
    float VL_F(Range x)
    float VG_F(Range x)
    float VT_F(Range x)
    float NP_F(Range x)
    float HC_F(Range x)
    float VC_F(Range x)
    float SD_F(Range x)
    
    Range DZ_INVF(float x)
    Range VR_INVF(float x)
    Range SW_INVF(float x)
    Range CZ_INVF(float x)
    Range ZT_INVF(float x)
    Range DR_INVF(float x)
    Range LR_INVF(float x)
    Range ZD_INVF(float x)
    Range DM_INVF(float x)
    Range RH_INVF(float x)
    Range PH_INVF(float x)
    Range XZ_INVF(float x)
    Range CD_INVF(float x)
    Range MZ_INVF(float x)
    Range MD_INVF(float x)
    Range ZE_INVF(float x)
    Range VE_INVF(float x)
    Range KD_INVF(float x)
    Range TI_INVF(float x)
    Range DX_INVF(float x)
    Range CH_INVF(float x)
    Range AH_INVF(float x)
    Range CV_INVF(float x)
    Range AV_INVF(float x)
    Range SQ_INVF(float x)
    Range VS_INVF(float x)
    Range VL_INVF(float x)
    Range VG_INVF(float x)
    Range VT_INVF(float x)
    Range NP_INVF(float x)
    Range HC_INVF(float x)
    Range VC_INVF(float x)
    Range SD_INVF(float x)
