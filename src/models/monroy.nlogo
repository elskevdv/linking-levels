;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                   Definition of Model Parameters                                   ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

globals [ Year Day half_per_patch speed_in_patches reference_T Boltz T Arrhenius current_incubation_period
  adult_size juvenile_size cocoon_size soil_moisture temperature
  
  saver_numAs                       ; list for saving the number of adults
  saver_numJs                       ; list for saving the number of juveniles
  saver_numCs                       ; list for saving the number of cocoons
   
  saver_massJs                      ; list for saving the number of juveniles
  saver_massAs                      ; list for saving the number of adults
  
  saved_nums                        ; list for saving the average number of model worms
  saved_mass                        ; list for saving the average total mass of model worms
  data_nums                         ; list for storing the average number of real worms
  ]

patches-own 
[ change_in_food                    ; change in food per timestep
  current_food                      ; food (g)
  func_response                     ; relative functional response
  ]

turtles-own 
[ mass                              ; individual mass (g)
  BMR                               ; energy cost of maintenance (kJ)
  ingestion_rate                    ; ingestion rate (g/day)
  
  energy_left                       ; intermediate daily energy reserve (kJ)
  energy_reserve_max                ; maximum energy reserve (kJ) (dependent on mass)
  energy_reserve                    ; amount of energy stored as tissue
  hatchlings                        ; number of hatchlings produced
  R                                 ; energy accumulated towards cocoon production (kJ)
  age
  embryonic_development
  ]

breed [ adults adult ]
breed [ juveniles juvenile ]
breed [ cocoons cocoon ]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                          Setup Functions                                           ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-interface
  set initial_number_cocoons 500
  set initial_number_juveniles 500
  set initial_number_adults 500
  set scape_size 1                  ; m
  set temperature 20                ; C
  set plot_year 2
end

to setup-basic-parameters
  set B_0 967                       ; kJ/day
  set activation_energy 0.25        ; eV
  set energy_tissue 7               ; kJ/g
  set energy_food 10.6              ; kJ/g
  set energy_synthesis 3.6          ; kJ/g

  set half_saturation_coef 3.5      ; g/0.01 m    
  set max_ingestion_rate 0.7        ; g/day 
  set mass_birth 0.011              ; g
  set mass_cocoon 0.015             ; g
  set mass_maximum  0.5             ; g
  
  set mass_sexual_maturity 0.25     ; g
  set growth_constant  0.177        ; g
  set max_reproduction_rate 0.182   ; kJ/g/day
  set speed 0.004                   ; m/day
  set incubation_period 23          ; day
end

to setup-ABC-parameters
  set B_0 1171.67                   ; kJ/day
  set activation_energy 0.237       ; eV
  set energy_tissue 5.366           ; kJ/g
  set energy_food 15.45             ; kJ/g
  set energy_synthesis 5.365        ; kJ/g

  set half_saturation_coef 5.534    ; g/0.01 m    
  set max_ingestion_rate 0.6859     ; g/day 
  set mass_birth 0.0124             ; g
  set mass_cocoon 0.0226            ; g
  set mass_maximum  0.5149          ; g
  
  set mass_sexual_maturity 0.1751   ; g
  set growth_constant  0.1668       ; g
  set max_reproduction_rate 0.0936  ; kJ/g/day
  set speed 0.008                   ; m/day
  set incubation_period 12.5        ; day
end


to setup
  clear-all
  
  set saved_nums (list [ ] [ ] [ ])
  set saved_mass (list [ ] [ ])
    
  let adult_nums (list 1680 1760 3810 432)
  let juv_nums (list 690 1790 3070 3992)
  let coc_nums (list 0 1260 3970 184)
  set data_nums (list adult_nums juv_nums coc_nums)
  
  set-default-shape adults "worm"
  set-default-shape juveniles "worm"
  set-default-shape cocoons "dot" 
  
  set adult_size 1
  set juvenile_size 0.8
  set cocoon_size 0.4
  
  set half_per_patch ((half_saturation_coef * scape_size) / 0.01) / count patches
  set speed_in_patches speed / (sqrt(scape_size) / sqrt(count patches))

  set reference_T 298.15            ; Kelvins
  set Boltz (8.62 * (10 ^ -5))      ; eV K-1
  
  set T 273.15 + temperature
  set Arrhenius (e ^ ((- activation_energy / Boltz ) * ((1 /  T ) - (1 / reference_T))))
  
  set Year 1
  set Day 1
  
  setup-patches
  setup-turtles
  
  reset-ticks
end

to setup-patches
  ask patches 
  [ update-patch ]
end

to setup-turtles
  create-adults initial_number_adults
  [ set color red
    set size adult_size
    setxy random-xcor random-ycor
    set mass mass_sexual_maturity
    set energy_reserve_max (mass / 2) * energy_tissue
    set energy_reserve energy_reserve_max ]
  
  create-juveniles initial_number_juveniles
  [ set color pink
    set size juvenile_size
    setxy random-xcor random-ycor
    set mass mass_birth
    set energy_reserve_max (mass / 2) * energy_tissue
    set energy_reserve energy_reserve_max ]
  
  create-cocoons initial_number_cocoons
  [ set color white
    set size cocoon_size
    setxy random-xcor random-ycor
    set mass mass_cocoon
    set energy_reserve_max (mass / 2) * energy_tissue
    set energy_reserve energy_reserve_max ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                            Run Functions                                           ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go
  if ticks = 731 [ stop ]
  start-saves
  update-environment
  calc-ingestion
  
  ask adults
  [ calc-mortality
    calc-assimilation
    calc-maintenance
    calc-reproduction
    calc-growth
    update-reserves
    move ]
  
  ask juveniles
  [ calc-mortality
    calc-assimilation
    calc-maintenance
    calc-growth
    transform-juvenile
    update-reserves
    move ]
  
  ask cocoons 
  [ set age age + 1
    calc-mortality
    calc-embryo-development
    calc-maintenance
    transform-cocoon ]
  
  ask patches
  [ update-patch ]
  
  save-nums
  plot-worms
  
  ifelse Day = 365 
  [ set Year Year + 1
    set Day 1 ]
  [ set Day Day + 1 ]
  
  tick
end

to go-to-end-of-data
  repeat (365 * plot_year) [ go ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                      Ingestion & Assimilation                                      ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-ingestion
  ask turtles with [ breed != cocoons ] [ set ingestion_rate (max_ingestion_rate * Arrhenius)
    * func_response * (mass ^ (2 / 3)) ]
  ask patches [ correct-ingestion-rate ]
end

to calc-assimilation
  set energy_left (ingestion_rate * energy_food)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                      Maintenance & Starvation                                      ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-maintenance
  set BMR B_0 * (mass ^ (3 / 4)) * e ^ (- activation_energy / (Boltz * T))     ; { Equation 1 }
  
  ifelse energy_left > BMR
  [ set energy_left energy_left - BMR ]
  [ set energy_reserve energy_reserve - (BMR - energy_left)
    set energy_left 0 ]
  
  if energy_reserve < (energy_reserve_max * 0.5) and breed != cocoons
  [ onset-starvation-strategy ]
end

to onset-starvation-strategy
  ifelse (energy_tissue + energy_synthesis > 0)
  [ set mass mass - (BMR / (energy_tissue + energy_synthesis)) ]
  [ die ]
  set energy_reserve energy_reserve + BMR
  
  if mass < mass_sexual_maturity 
  [ set breed juveniles
    set color pink
    set size juvenile_size ]
  
  if mass < mass_birth 
  [ die ]  
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                            Reproduction                                            ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-reproduction
  let max_R (max_reproduction_rate * Arrhenius) * mass
  
  let energy_for_R min list energy_left max_R
  set energy_left energy_left - energy_for_R
  
  if (energy_for_R < max_R and energy_reserve > 0.5 * energy_reserve_max)
  [ set energy_reserve energy_reserve - (max_R - energy_for_R)
    set energy_for_R max_R ] 
  
  set R R + energy_for_R
  
  if R >= (mass_cocoon * (energy_tissue + energy_synthesis))
  [ reproduce ]
end

to reproduce
  hatch-cocoons 1
  [ set color white
    set size cocoon_size
    set mass mass_cocoon
    set ingestion_rate 0
    set energy_left 0
    set energy_reserve_max (mass_cocoon * energy_tissue)
    set energy_reserve energy_reserve_max
    set hatchlings 0
    set R 0
    set age 0
    set embryonic_development 0 ]
  
  set R (R - (mass_cocoon * (energy_tissue + energy_synthesis))) 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                                Growth                                              ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-growth
  let max_G (growth_constant * Arrhenius) * (mass_maximum ^ (1 / 3) * mass ^ (2 / 3) - mass)
  let energy_for_G min list (max_G * (energy_tissue + energy_synthesis)) energy_left
  
  let to_grow 0
  
  if (max_G > 0 and (energy_tissue + energy_synthesis) > 0)
  [ set to_grow (energy_for_G / (max_G * (energy_tissue + energy_synthesis))) * max_G ]  
  
  if (mass + to_grow) < mass_maximum
  [ set mass mass + to_grow
    set energy_left energy_left - energy_for_G ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                           Energy Reserves                                          ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-reserves
  set energy_reserve_max (mass / 2) * energy_tissue
  
  if energy_left > 0
  [ ifelse (energy_tissue + energy_synthesis > 0)
    [ set energy_reserve energy_reserve + (energy_left *
        (energy_tissue / (energy_tissue + energy_synthesis))) ]
    [ die ] ]
  
  ;;; a maximum threshold of 50% mass * energy content of flesh is set for energy reserves   
  if energy_reserve > energy_reserve_max and breed != cocoons 
  [ set energy_reserve energy_reserve_max ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                              Movement                                              ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to move   
    rt random 90
    lt random 90
    fd speed_in_patches
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                       Life Stage Development                                       ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-mortality
  let mort (12.7 - 0.001 * soil_moisture - 0.0861 * T + 0.000009 *
    (soil_moisture ^ 2) + 0.000147 * (T ^ 2)) / 10
  if random-float 1 < mort
  [ die ]
end

to calc-embryo-development
   set embryonic_development embryonic_development + ((1 / current_incubation_period) * 100)
end

to transform-cocoon
  if age >= incubation_period and embryonic_development >= 100
  [ set breed juveniles
    set color pink
    set size juvenile_size ] 
end

to transform-juvenile
  if mass >= mass_sexual_maturity 
  [ set breed adults
    set color red
    set size adult_size ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                           Patch Functions                                          ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-patch
  if current_food > 0 
  [ set current_food max list (current_food - sum [ ingestion_rate ]
      of turtles-here with [ breed != cocoons ]) 0
    set func_response current_food / (current_food + half_per_patch) ]

  ifelse current_food > 0 
  [ set pcolor scale-color green current_food (150 * 2) (0 - 150 * 0.5) ]
  [ set pcolor brown ]
end

to correct-ingestion-rate
  if current_food = 0
  [ ask turtles-here with [ breed != cocoons ] [ set ingestion_rate 0 ] ]
  
  if sum [ ingestion_rate ] of turtles-here with [ breed != cocoons ] > current_food
  [ ask turtles-here with [ breed != cocoons ] [ set ingestion_rate current_food /
      count turtles-here with [ breed != cocoons ] ]
    set current_food 0 ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                        Experiment Simulation                                       ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-environment
  if Day >= 1 and Day < 92 
  [ set temperature random-normal 9 0.9
    set soil_moisture random-normal 60 6 ]
  
  if Day >= 92 and Day < 183
  [ set temperature random-normal 8 0.8
    set soil_moisture random-normal 40 4 ]
  
  if Day >= 183 and Day < 274
  [ set temperature random-normal 17 1.7
    set soil_moisture random-normal 70 7 ]
  
  if Day >= 274 and Day < 365
  [ set temperature random-normal 21 2.1
    set soil_moisture random-normal 45 4.5 ]
  
  set T 273.15 + temperature
  set Arrhenius (e ^ ((- activation_energy / Boltz ) * ((1 /  T ) - (1 / reference_T))))
  set current_incubation_period incubation_period * e ^ ((- activation_energy / Boltz ) *
    ((1 /  reference_T) - (1 / T)))
  
  if Day = 1
  [ ask patches [ set current_food abs (random-normal (100 / count patches) (10 / count patches)) * 100 ] ]
  
  if Day = 92
  [ ask patches [ set current_food abs (random-normal (150 / count patches) (15 / count patches)) * 100 ] ]
  
  if Day = 183
  [ ask patches [ set current_food abs (random-normal (500 / count patches) (50 / count patches)) * 100 ] ]

  if Day = 274
  [ ask patches [ set current_food abs (random-normal (50 / count patches) (5 / count patches)) * 100 ] ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                             Reporters                                              ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report total-adult-mass
  ifelse any? adults
  [ report sum [ mass ] of adults ]
  [ report 0 ]
end

to-report total-juvenile-mass
  ifelse any? juveniles
  [ report sum [ mass ] of juveniles ]
  [ report 0 ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                           Plot Functions                                           ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to start-saves
  if Year = plot_year and (Day = 1 or Day = 91 or Day = 181 or Day = 271)
  [ set saver_numAs [ ] 
    set saver_numJs [ ]
    set saver_numCs [ ]

    set saver_massAs [ ]   
    set saver_massJs [ ] ]

end

to save-nums
  if Year = plot_year
  [ set saver_numAs lput count adults saver_numAs
    set saver_numJs lput count juveniles saver_numJs
    set saver_numCs lput count cocoons saver_numCs
    
    set saver_massAs lput total-adult-mass saver_massAs
    set saver_massJs lput total-juvenile-mass saver_massJs

    if Day = 90 or Day = 180 or Day = 270 or Day = 364
    [ set saved_nums replace-item 0 saved_nums lput mean saver_numAs item 0 saved_nums
      set saved_nums replace-item 1 saved_nums lput mean saver_numJs item 1 saved_nums
      set saved_nums replace-item 2 saved_nums lput mean saver_numCs item 2 saved_nums
      
      set saved_mass replace-item 0 saved_mass lput mean saver_massAs item 0 saved_mass
      set saved_mass replace-item 1 saved_mass lput mean saver_massJs item 1 saved_mass ] ]
end

to plot-worms
  if Year = plot_year
  [ set-current-plot "# of Adults"
    set-current-plot-pen "Model_D"
    plot count adults
    set-current-plot-pen "Model_S"
    if Day = 90 [ plotxy  (0 +  90 - 46) item 0 item 0 saved_nums ]
    if Day = 180 [ plotxy (0 + 180 - 46) item 1 item 0 saved_nums ]
    if Day = 270 [ plotxy (0 + 270 - 46) item 2 item 0 saved_nums ]
    if Day = 364 [ plotxy (0 + 364 - 46) item 3 item 0 saved_nums ]
    
    set-current-plot-pen "Data"
    if Day = 90 [ plotxy  (0 +  90 - 46) item 0 item 0 data_nums ]
    if Day = 180 [ plotxy (0 + 180 - 46) item 1 item 0 data_nums ]
    if Day = 270 [ plotxy (0 + 270 - 46) item 2 item 0 data_nums ]
    if Day = 364 [ plotxy (0 + 364 - 46) item 3 item 0 data_nums ]
    
    set-current-plot "# of Juveniles"
    set-current-plot-pen "Model_D"
    plot count juveniles
    set-current-plot-pen "Model_S"
    if Day = 90 [ plotxy  (0 +  90 - 46) item 0 item 1 saved_nums ]
    if Day = 180 [ plotxy (0 + 180 - 46) item 1 item 1 saved_nums ]
    if Day = 270 [ plotxy (0 + 270 - 46) item 2 item 1 saved_nums ]
    if Day = 364 [ plotxy (0 + 364 - 46) item 3 item 1 saved_nums ]
    
    set-current-plot-pen "Data"
    if Day = 90 [ plotxy  (0 +  90 - 46) item 0 item 1 data_nums ]
    if Day = 180 [ plotxy (0 + 180 - 46) item 1 item 1 data_nums ]
    if Day = 270 [ plotxy (0 + 270 - 46) item 2 item 1 data_nums ]
    if Day = 364 [ plotxy (0 + 364 - 46) item 3 item 1 data_nums ]
    
    set-current-plot "# of Cocoons"
    set-current-plot-pen "Model_D"
    plot count cocoons
    set-current-plot-pen "Model_S"
    if Day = 90 [ plotxy  (0 +  90 - 46) item 0 item 2 saved_nums ]
    if Day = 180 [ plotxy (0 + 180 - 46) item 1 item 2 saved_nums ]
    if Day = 270 [ plotxy (0 + 270 - 46) item 2 item 2 saved_nums ]
    if Day = 364 [ plotxy (0 + 364 - 46) item 3 item 2 saved_nums ]
    
    set-current-plot-pen "Data"
    if Day = 90 [ plotxy  (0 +  90 - 46) item 0 item 2 data_nums ]
    if Day = 180 [ plotxy (0 + 180 - 46) item 1 item 2 data_nums ]
    if Day = 270 [ plotxy (0 + 270 - 46) item 2 item 2 data_nums ]
    if Day = 364 [ plotxy (0 + 364 - 46) item 3 item 2 data_nums ] ]
end
@#$#@#$#@
GRAPHICS-WINDOW
263
22
593
373
-1
-1
20.0
1
10
1
1
1
0
1
1
1
0
15
0
15
1
1
1
ticks
30.0

BUTTON
990
271
1054
304
NIL
Setup\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1129
271
1192
304
NIL
Go\n
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
32
64
246
97
initial_number_juveniles
initial_number_juveniles
0
5000
500
10
1
NIL
HORIZONTAL

PLOT
262
381
594
524
Total Food
time (ticks)
food (g)
0.0
80.0
0.0
200.0
true
true
"" ""
PENS
"total food" 1.0 0 -16777216 true "" "plot sum [ current_food ] of patches"
"total ingestion" 1.0 0 -7500403 true "" "plot sum [ ingestion_rate ] of turtles with [ breed != cocoons ]"

SLIDER
31
105
247
138
initial_number_adults
initial_number_adults
0
5000
500
10
1
NIL
HORIZONTAL

PLOT
604
464
981
645
# of Cocoons
time
# of worms
0.0
2.0
0.0
0.05
true
false
"" ""
PENS
"Model_D" 1.0 1 -11221820 true "" ""
"Model_S" 1.0 0 -9276814 true "" ""
"Data" 1.0 0 -16777216 true "" ""

BUTTON
1060
271
1123
304
Go
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
991
61
1061
121
B_0
1171.67
1
0
Number

INPUTBOX
1067
62
1177
122
activation_energy
0.237
1
0
Number

INPUTBOX
1184
62
1268
122
energy_tissue
5.366
1
0
Number

INPUTBOX
1391
62
1493
122
energy_synthesis
5.365
1
0
Number

INPUTBOX
1118
130
1237
190
max_ingestion_rate
0.6859
1
0
Number

INPUTBOX
992
130
1111
190
half_saturation_coef
5.534
1
0
Number

INPUTBOX
1245
131
1317
191
mass_birth
0.0124
1
0
Number

INPUTBOX
991
200
1121
260
mass_sexual_maturity
0.1751
1
0
Number

INPUTBOX
1411
133
1513
193
mass_maximum
0.5149
1
0
Number

INPUTBOX
1323
132
1405
192
mass_cocoon
0.0226
1
0
Number

INPUTBOX
1129
200
1228
260
growth_constant
0.1668
1
0
Number

INPUTBOX
1237
199
1371
259
max_reproduction_rate
0.0936
1
0
Number

BUTTON
1119
311
1253
344
Set Basic Parameters
setup-basic-parameters
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
990
311
1111
344
Set Basic Interface
setup-interface
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1276
62
1385
122
energy_food
15.45
1
0
Number

INPUTBOX
1379
200
1451
260
speed
0.0080
1
0
Number

INPUTBOX
30
145
129
205
scape_size
1
1
0
Number

SLIDER
33
24
246
57
initial_number_cocoons
initial_number_cocoons
0
5000
500
10
1
NIL
HORIZONTAL

INPUTBOX
1459
201
1567
261
incubation_period
12.5
1
0
Number

PLOT
263
536
594
686
Daily Mortality Rates
NIL
NIL
0.0
10.0
0.0
0.01
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot (12.7 - 0.001 * soil_moisture - 0.0861 * T + 0.000009 * (soil_moisture ^ 2) + 0.000147 * (T ^ 2)) / 10"

PLOT
605
270
978
452
# of Juveniles
time
# of worms
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"Model_D" 1.0 1 -2064490 true "" ""
"Model_S" 1.0 0 -9276814 true "" ""
"Data" 1.0 0 -16777216 true "" ""

PLOT
604
77
980
259
# of Adults
time
# of worms
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"Model_D" 1.0 1 -2674135 true "" ""
"Model_S" 1.0 0 -9276814 true "" ""
"Data" 1.0 0 -16777216 true "" ""

INPUTBOX
136
145
247
205
plot_year
2
1
0
Number

BUTTON
1381
311
1520
344
Go to End of Data
go-to-end-of-data
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
135
213
246
258
NIL
temperature
2
1
11

MONITOR
30
213
129
258
NIL
soil_moisture
2
1
11

MONITOR
604
23
661
68
NIL
Year
17
1
11

MONITOR
670
23
727
68
NIL
Day
17
1
11

BUTTON
1260
311
1374
344
Set ABC Parameters
setup-ABC-parameters
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
608
656
978
691
Note: Plots will start appearing once the simulation reaches \"plot_year\" as set at left.
11
0.0
1

@#$#@#$#@
## CREDITS AND REFERENCES

Copyright (C) 2015 Elske van der Vaart, Alice Johnston, Richard Sibly
<elskevdv at gmail.com>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

This program accompanies the following paper:
van der Vaart, E., Johnston, A.S.A., & Sibly, R.M. "Predicting how many animals will be where: How to build, calibrate and evaluate individual-based models"
(2015) Ecological Modelling.

The earthworm IBM it runs was originally described here:
Johnston, A.S.A., Hodson, M.E., Thorbek, P., Alvarez, T. & Sibly, R.M. "An energy budget agent-based model of earthworm populations and its application to study the effects of pesticides" (2014) Ecological Modelling, 280, 5 - 17.

The empirical data it plots was taken from:
Monroy, F., Aira, M., Domínguez, J., & Velando, A. "Seasonal population dynamics of Eisenia fetida (Savigny, 1826) (Oligochaeta, Lumbricidae) in the field" (2006)
Comptes Rendus Biologies, 329, 912 - 915.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

worm
true
0
Polygon -7500403 true true 165 210 165 225 135 255 105 270 90 270 75 255 75 240 90 210 120 195 135 165 165 135 165 105 150 75 150 60 135 60 120 45 120 30 135 15 150 15 180 30 180 45 195 45 210 60 225 105 225 135 210 150 210 165 195 195 180 210

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Temp16.5" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>ticks = 365</exitCondition>
    <metric>count juveniles + count adults + count cocoons</metric>
    <metric>count adults</metric>
    <metric>sum [hatchlings] of turtles with [breed = adults]</metric>
    <metric>count juveniles</metric>
    <metric>sum [mass] of turtles with [breed = juveniles]</metric>
  </experiment>
  <experiment name="getal:random" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="189"/>
    <exitCondition>ticks = 189</exitCondition>
    <metric>sum [mass] of turtles with [breed = juveniles]</metric>
    <metric>sum [mass] of turtles with [breed = adults]</metric>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="180"/>
    <exitCondition>ticks = 180</exitCondition>
    <metric>sum [hatchlings] of turtles with [breed = adults]</metric>
  </experiment>
  <experiment name="experiment" repetitions="25" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>ticks = 189</exitCondition>
    <metric>sum [mass] of turtles</metric>
    <enumeratedValueSet variable="Initial_number_cocoons">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial_number_adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food_dynamics">
      <value value="&quot;depleting&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="soil_moisture">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food_density_patch">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial_number_juveniles">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Temperature">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temp?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="moisture?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Elske's Experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup
repeat 3650 [ go ]</setup>
    <go>go</go>
    <metric>count cocoons</metric>
    <metric>count juveniles</metric>
    <metric>count adults</metric>
    <enumeratedValueSet variable="B_0">
      <value value="967"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scape_size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="speed">
      <value value="0.0040"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial_number_adults">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mass_sexual_maturity">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incubation_period">
      <value value="23"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="energy_intake">
      <value value="10.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mass_cocoon">
      <value value="0.015"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_ingestion_rate">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mass_birth">
      <value value="0.011"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial_number_cocoons">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="energy_flesh">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="half_saturation_coef">
      <value value="3.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_reproduction_rate">
      <value value="0.182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="energy_synthesis">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="soil_moisture">
      <value value="58.14355285443175"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="growth_constant">
      <value value="0.177"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="8.485796762071958"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="activation_energy">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mass_maximum">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial_number_juveniles">
      <value value="500"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
