	  $%  k   k820309              13.0        �q4S                                                                                                           
       /Users/adamwolf/Research/LM3-PPA-r44/src/land_lad2/shared/land_debug.F90 LAND_DEBUG_MOD              LAND_DEBUG_INIT LAND_DEBUG_END SET_CURRENT_POINT GET_CURRENT_POINT CURRENT_I CURRENT_J CURRENT_K CURRENT_FACE IS_WATCH_POINT GET_WATCH_POINT CHECK_TEMP_RANGE CHECK_CONSERVATION WATER_CONS_TOL CARBON_CONS_TOL gen@DPRI                      @                              
       LOGUNIT GET_UNIT ERROR_MESG FATAL CHECK_NML_ERROR                      @                              
       TIME_TYPE GET_DATE                                                        u #DEBUG_PRINTOUT_R0D    #DEBUG_PRINTOUT_I0D    #DEBUG_PRINTOUT_L0D    #DEBUG_PRINTOUT_R1D    #DEBUG_PRINTOUT_I1D    #DEBUG_PRINTOUT_R2D    #         @     @X                                                #DEBUG_PRINTOUT_R0D%LEN_TRIM    #DEBUG_PRINTOUT_R0D%TRIM    #DESCRIPTION    #VALUE                  @                                 LEN_TRIM               @                                 TRIM           
  @@                                                 1           
   @                                   
      #         @     @X                                                #DEBUG_PRINTOUT_I0D%LEN_TRIM 	   #DEBUG_PRINTOUT_I0D%TRIM 
   #DESCRIPTION    #VALUE                  @                            	     LEN_TRIM               @                            
     TRIM           
  @@                                                 1           
   @                                         #         @     @X                                                #DEBUG_PRINTOUT_L0D%LEN_TRIM    #DEBUG_PRINTOUT_L0D%TRIM    #DESCRIPTION    #VALUE                  @                                 LEN_TRIM               @                                 TRIM           
  @@                                                 1           
   @                                         #         @     @X                                                #DEBUG_PRINTOUT_R1D%SIZE    #DEBUG_PRINTOUT_R1D%LEN_TRIM    #DEBUG_PRINTOUT_R1D%TRIM    #DESCRIPTION    #VALUES                  @                                 SIZE               @                                 LEN_TRIM               @                                 TRIM           
  @@                                                 1           
 @@                                                 
              &                                           #         @     @X                                                #DEBUG_PRINTOUT_I1D%LEN_TRIM    #DEBUG_PRINTOUT_I1D%TRIM    #DESCRIPTION    #VALUES                  @                                 LEN_TRIM               @                                 TRIM           
  @@                                                 1           
   @                                                               &                                           #         @     @X                                                #DEBUG_PRINTOUT_R2D%LEN_TRIM    #DEBUG_PRINTOUT_R2D%TRIM    #DESCRIPTION     #VALUES !                 @                                 LEN_TRIM               @                                 TRIM           
  @@                                                  1           
   @                              !                   
              &                   &                                                          �  !@                           "     '                    #SECONDS #   #DAYS $   #TICKS %   #DUMMY &                � D                              #                                � D                              $                               � D                              %                               � D                              &                                   !                            '                                       
               10%         @     !                           (                            #         @       !                           )                   #ERROR_MESG%TRIM *   #ROUTINE +   #MESSAGE ,   #LEVEL -                 @                            *     TRIM           
   @                             +                    1           
   @                             ,                    1           
   @                              -                            !                            .                                                      2%         @      !                          /                          #CHECK_NML_ERROR%TRIM 0   #IOSTAT 1   #NML_NAME 2                 @                            0     TRIM           
   @                              1                     
   @                             2                    1 #         @       !                           3                	   #GET_DATE%PRESENT 4   #TIME 5   #YEAR 6   #MONTH 7   #DAY 8   #HOUR 9   #MINUTE :   #SECOND ;   #TICK <   #ERR_MSG =                 @                            4     PRESENT           
   @                              5                   #TIME_TYPE "               @                              6                        @                              7                        @                              8                        @                              9                        @                              :                        @                              ;                        @                              <                        @                             =                     1 #         @                                   >                    #LAND_DEBUG_INIT%TRIM ?                 @                            ?     TRIM #         @                                   @                     #         @                                   A                    #I B   #J C   #K D             
   @                              B                     
   @                              C                     
   @                              D           #         @                                   E                   #GET_CURRENT_POINT%PRESENT F   #I G   #J H   #K I   #FACE J                 @                            F     PRESENT           F @@                              G                      F @@                              H                      F @@                              I                      F @@                              J            %         @                                K                            %         @                                L                            %         @                                M                            %         @                                N                            %         @                                O                            #         @                                   P                   #GET_WATCH_POINT%PRESENT Q   #I R   #J S   #K T   #FACE U                 @                            Q     PRESENT           F @@                              R                      F @@                              S                      F @@                              T                      F @@                              U            #         @                                   V                   #CHECK_TEMP_RANGE%TRIM W   #TEMP X   #TAG Y   #VARNAME Z   #TIME [                 @                            W     TRIM           
   @                              X     
                
  @@                             Y                    1           
  @@                             Z                    1           
  @@                              [                   #TIME_TYPE "   #         @                                   \                   #CHECK_CONSERVATION%ABS ]   #CHECK_CONSERVATION%PRESENT ^   #CHECK_CONSERVATION%TRIM _   #TAG `   #SUBSTANCE a   #D1 b   #D2 c   #TOLERANCE d   #TIME e   #SEVERITY f                 @                            ]     ABS               @                            ^     PRESENT               @                            _     TRIM           
  @@                             `                    1           
  @@                             a                    1           
   @                              b     
                
   @                              c     
                
   @                              d     
                
  @@                              e                   #TIME_TYPE "             
 @@                              f                     D@                                g     
                 D@                                h     
          �   `      fn#fn $      �   b   uapp(LAND_DEBUG_MOD    �  r   J  UTILITIES_MOD !   [  S   J  TIME_MANAGER_MOD    �  �       gen@DPRI #   ~  �      DEBUG_PRINTOUT_R0D ,      A      DEBUG_PRINTOUT_R0D%LEN_TRIM (   a  =      DEBUG_PRINTOUT_R0D%TRIM /   �  L   a   DEBUG_PRINTOUT_R0D%DESCRIPTION )   �  @   a   DEBUG_PRINTOUT_R0D%VALUE #   *  �      DEBUG_PRINTOUT_I0D ,   �  A      DEBUG_PRINTOUT_I0D%LEN_TRIM (     =      DEBUG_PRINTOUT_I0D%TRIM /   J  L   a   DEBUG_PRINTOUT_I0D%DESCRIPTION )   �  @   a   DEBUG_PRINTOUT_I0D%VALUE #   �  �      DEBUG_PRINTOUT_L0D ,   x  A      DEBUG_PRINTOUT_L0D%LEN_TRIM (   �  =      DEBUG_PRINTOUT_L0D%TRIM /   �  L   a   DEBUG_PRINTOUT_L0D%DESCRIPTION )   B  @   a   DEBUG_PRINTOUT_L0D%VALUE #   �  �      DEBUG_PRINTOUT_R1D (   B	  =      DEBUG_PRINTOUT_R1D%SIZE ,   	  A      DEBUG_PRINTOUT_R1D%LEN_TRIM (   �	  =      DEBUG_PRINTOUT_R1D%TRIM /   �	  L   a   DEBUG_PRINTOUT_R1D%DESCRIPTION *   I
  �   a   DEBUG_PRINTOUT_R1D%VALUES #   �
  �      DEBUG_PRINTOUT_I1D ,   x  A      DEBUG_PRINTOUT_I1D%LEN_TRIM (   �  =      DEBUG_PRINTOUT_I1D%TRIM /   �  L   a   DEBUG_PRINTOUT_I1D%DESCRIPTION *   B  �   a   DEBUG_PRINTOUT_I1D%VALUES #   �  �      DEBUG_PRINTOUT_R2D ,   q  A      DEBUG_PRINTOUT_R2D%LEN_TRIM (   �  =      DEBUG_PRINTOUT_R2D%TRIM /   �  L   a   DEBUG_PRINTOUT_R2D%DESCRIPTION *   ;  �   a   DEBUG_PRINTOUT_R2D%VALUES +   �  }       TIME_TYPE+TIME_MANAGER_MOD ;   \  H   %   TIME_TYPE%SECONDS+TIME_MANAGER_MOD=SECONDS 5   �  H   %   TIME_TYPE%DAYS+TIME_MANAGER_MOD=DAYS 7   �  H   %   TIME_TYPE%TICKS+TIME_MANAGER_MOD=TICKS 7   4  H   %   TIME_TYPE%DUMMY+TIME_MANAGER_MOD=DUMMY &   |  r       LOGUNIT+UTILITIES_MOD '   �  P       GET_UNIT+UTILITIES_MOD )   >  �       ERROR_MESG+UTILITIES_MOD 3   �  =      ERROR_MESG%TRIM+UTILITIES_MOD=TRIM 1   �  L   e   ERROR_MESG%ROUTINE+UTILITIES_MOD 1   I  L   e   ERROR_MESG%MESSAGE+UTILITIES_MOD /   �  @   e   ERROR_MESG%LEVEL+UTILITIES_MOD $   �  q       FATAL+UTILITIES_MOD .   F  �       CHECK_NML_ERROR+UTILITIES_MOD 8   �  =      CHECK_NML_ERROR%TRIM+UTILITIES_MOD=TRIM 5     @   e   CHECK_NML_ERROR%IOSTAT+UTILITIES_MOD 7   G  L   e   CHECK_NML_ERROR%NML_NAME+UTILITIES_MOD *   �  �       GET_DATE+TIME_MANAGER_MOD :   R  @      GET_DATE%PRESENT+TIME_MANAGER_MOD=PRESENT /   �  W   e   GET_DATE%TIME+TIME_MANAGER_MOD /   �  @   e   GET_DATE%YEAR+TIME_MANAGER_MOD 0   )  @   e   GET_DATE%MONTH+TIME_MANAGER_MOD .   i  @   e   GET_DATE%DAY+TIME_MANAGER_MOD /   �  @   e   GET_DATE%HOUR+TIME_MANAGER_MOD 1   �  @   e   GET_DATE%MINUTE+TIME_MANAGER_MOD 1   )  @   e   GET_DATE%SECOND+TIME_MANAGER_MOD /   i  @   e   GET_DATE%TICK+TIME_MANAGER_MOD 2   �  L   e   GET_DATE%ERR_MSG+TIME_MANAGER_MOD     �  b       LAND_DEBUG_INIT %   W  =      LAND_DEBUG_INIT%TRIM    �  H       LAND_DEBUG_END "   �  ]       SET_CURRENT_POINT $   9  @   a   SET_CURRENT_POINT%I $   y  @   a   SET_CURRENT_POINT%J $   �  @   a   SET_CURRENT_POINT%K "   �  �       GET_CURRENT_POINT *     @      GET_CURRENT_POINT%PRESENT $   �  @   a   GET_CURRENT_POINT%I $   �  @   a   GET_CURRENT_POINT%J $   ?  @   a   GET_CURRENT_POINT%K '     @   a   GET_CURRENT_POINT%FACE    �  P       CURRENT_I      P       CURRENT_J    _  P       CURRENT_K    �  P       CURRENT_FACE    �  P       IS_WATCH_POINT     O  �       GET_WATCH_POINT (   �  @      GET_WATCH_POINT%PRESENT "     @   a   GET_WATCH_POINT%I "   S  @   a   GET_WATCH_POINT%J "   �  @   a   GET_WATCH_POINT%K %   �  @   a   GET_WATCH_POINT%FACE !     �       CHECK_TEMP_RANGE &   �  =      CHECK_TEMP_RANGE%TRIM &   �  @   a   CHECK_TEMP_RANGE%TEMP %      L   a   CHECK_TEMP_RANGE%TAG )   i   L   a   CHECK_TEMP_RANGE%VARNAME &   �   W   a   CHECK_TEMP_RANGE%TIME #   !  �       CHECK_CONSERVATION '   �!  <      CHECK_CONSERVATION%ABS +   8"  @      CHECK_CONSERVATION%PRESENT (   x"  =      CHECK_CONSERVATION%TRIM '   �"  L   a   CHECK_CONSERVATION%TAG -   #  L   a   CHECK_CONSERVATION%SUBSTANCE &   M#  @   a   CHECK_CONSERVATION%D1 &   �#  @   a   CHECK_CONSERVATION%D2 -   �#  @   a   CHECK_CONSERVATION%TOLERANCE (   $  W   a   CHECK_CONSERVATION%TIME ,   d$  @   a   CHECK_CONSERVATION%SEVERITY    �$  @       WATER_CONS_TOL     �$  @       CARBON_CONS_TOL 