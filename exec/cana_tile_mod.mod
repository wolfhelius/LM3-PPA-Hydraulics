	    6   k820309              13.0        !�R                                                                                                           
       /Users/adamwolf/Research/LM3-PPA-r44/src/land_lad2/canopy_air/cana_tile.F90 CANA_TILE_MOD              CANA_PROG_TYPE CANA_TILE_TYPE DELETE_CANA_TILE CANA_TILES_CAN_BE_MERGED MERGE_CANA_TILES GET_CANA_TILE_TAG CANA_IS_SELECTED CANA_TILE_STOCK_PE CANA_TILE_CARBON CANA_TILE_HEAT CANOPY_AIR_MASS CANOPY_AIR_MASS_FOR_TRACERS CPW gen@NEW_CANA_TILE                      @                              
       TILE_SELECTOR_TYPE                      @                              
       CP_AIR TFREEZE                      @                              
       MOL_C MOL_CO2                                                        u #CANA_TILE_CTOR    #CANA_TILE_COPY_CTOR    &         @   @X                                                        #CANA_TILE_TYPE    &         @   @X                                                       #CANA    #CANA_TILE_TYPE              
   @                                                 #CANA_TILE_TYPE                      !@                               '�                    #NAME 	   #LONG_NAME 
   #MASK    #TAG    #IDATA1    #IDATA2    #RDATA1    #RDATA2                � $                             	                                                                                          �              C                                                � $                             
     �                                                                             �       �              C                                                                                                                                                               �$                                         �                            &                                                                                  y                                                           � $                                   �                                                                                           0                � $                                   �                                                                                           0                � $                                   �                                                                                           0                � $                                   �                                                                                           0                � $                                   �                                                                                           0                     !                                 
                   
                  ��Q�e�@                         !                                 
                 
                 ��(\�q@        273.16                 !                                 
                 
                 �~j�t��?        12.0E-3                 !                                 
                 
                 E���x��?        44.00995E-3                  @                               '                    #T    #Q    #CO2                 � $                                              
                � $                                             
                � $                                             
                     @                                '                    #PROG                 � $                                                         #CANA_PROG_TYPE    #         @                                                       #CANA              D @                                                  #CANA_TILE_TYPE    %         @                                                            #CANA1    #CANA2              
   @                                                 #CANA_TILE_TYPE              
   @                                                 #CANA_TILE_TYPE    #         @                                                       #CANA1     #W1 !   #CANA2 "   #W2 #             
   @                                                  #CANA_TILE_TYPE              
   @                              !     
                
D  @                              "                    #CANA_TILE_TYPE              
   @                              #     
      %         @                                 $                           #CANA %             
   @                              %                   #CANA_TILE_TYPE    %         @                                 &                           #CANA '   #SEL (             
   @                              '                   #CANA_TILE_TYPE              
   @                              (     �              #TILE_SELECTOR_TYPE    #         @                                   )                    #CANA *   #TWD_LIQ +   #TWD_SOL ,             
   @                              *                   #CANA_TILE_TYPE              D  @                              +     
                 D  @                              ,     
       %         @                                 -                    
       #CANA .             
   @                              .                   #CANA_TILE_TYPE    %         @                                 /                    
       #CANA 0             
   @                              0                   #CANA_TILE_TYPE                                               1     
                                                  2     
                                                  3     
          �   b      fn#fn #       b   uapp(CANA_TILE_MOD (     S   J  LAND_TILE_SELECTORS_MOD    V  O   J  CONSTANTS_MOD #   �  N   J  LAND_CONSTANTS_MOD "   �  m       gen@NEW_CANA_TILE    `  d      CANA_TILE_CTOR $   �  n      CANA_TILE_COPY_CTOR )   2  \   a   CANA_TILE_COPY_CTOR%CANA ;   �  �       TILE_SELECTOR_TYPE+LAND_TILE_SELECTORS_MOD @   :  �   a   TILE_SELECTOR_TYPE%NAME+LAND_TILE_SELECTORS_MOD E     =  a   TILE_SELECTOR_TYPE%LONG_NAME+LAND_TILE_SELECTORS_MOD @   D  �   a   TILE_SELECTOR_TYPE%MASK+LAND_TILE_SELECTORS_MOD ?   8  �   a   TILE_SELECTOR_TYPE%TAG+LAND_TILE_SELECTORS_MOD B   �  �   a   TILE_SELECTOR_TYPE%IDATA1+LAND_TILE_SELECTORS_MOD B   �	  �   a   TILE_SELECTOR_TYPE%IDATA2+LAND_TILE_SELECTORS_MOD B   '
  �   a   TILE_SELECTOR_TYPE%RDATA1+LAND_TILE_SELECTORS_MOD B   �
  �   a   TILE_SELECTOR_TYPE%RDATA2+LAND_TILE_SELECTORS_MOD %   q  p       CP_AIR+CONSTANTS_MOD &   �  v       TFREEZE+CONSTANTS_MOD )   W  w       MOL_C+LAND_CONSTANTS_MOD +   �  {       MOL_CO2+LAND_CONSTANTS_MOD    I  g       CANA_PROG_TYPE !   �  H   a   CANA_PROG_TYPE%T !   �  H   a   CANA_PROG_TYPE%Q #   @  H   a   CANA_PROG_TYPE%CO2    �  Z       CANA_TILE_TYPE $   �  d   a   CANA_TILE_TYPE%PROG !   F  R       DELETE_CANA_TILE &   �  \   a   DELETE_CANA_TILE%CANA )   �  f       CANA_TILES_CAN_BE_MERGED /   Z  \   a   CANA_TILES_CAN_BE_MERGED%CANA1 /   �  \   a   CANA_TILES_CAN_BE_MERGED%CANA2 !     n       MERGE_CANA_TILES '   �  \   a   MERGE_CANA_TILES%CANA1 $   �  @   a   MERGE_CANA_TILES%W1 '     \   a   MERGE_CANA_TILES%CANA2 $   x  @   a   MERGE_CANA_TILES%W2 "   �  Z       GET_CANA_TILE_TAG '     \   a   GET_CANA_TILE_TAG%CANA !   n  c       CANA_IS_SELECTED &   �  \   a   CANA_IS_SELECTED%CANA %   -  `   a   CANA_IS_SELECTED%SEL #   �  l       CANA_TILE_STOCK_PE (   �  \   a   CANA_TILE_STOCK_PE%CANA +   U  @   a   CANA_TILE_STOCK_PE%TWD_LIQ +   �  @   a   CANA_TILE_STOCK_PE%TWD_SOL !   �  Z       CANA_TILE_CARBON &   /  \   a   CANA_TILE_CARBON%CANA    �  Z       CANA_TILE_HEAT $   �  \   a   CANA_TILE_HEAT%CANA     A  @       CANOPY_AIR_MASS ,   �  @       CANOPY_AIR_MASS_FOR_TRACERS    �  @       CPW 