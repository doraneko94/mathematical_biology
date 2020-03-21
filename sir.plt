
set termoption dashed
set termoption enhanced
set tics front
reset
set xlabel "Year"
set ylabel "S, I, R"
set title "SIR model"
set key at  graph 1.0000000000000000e0, graph 6.9999999999999996e-1 right top vertical noreverse noinvert
unset xzeroaxis 
unset logscale x
set xdata
set mxtics 1

set xrange [0.000000000000e0:2.500000000000e1]
set xtics autofreq scale 5.000000000000e-1,5.000000000000e-1
unset yzeroaxis 
unset logscale y
set ydata
set mytics 1

set yrange [0.000000000000e0:1.000000000000e2]
set ytics autofreq scale 5.000000000000e-1,5.000000000000e-1

unset logscale cb
set cbdata
set mcbtics 1

set cbrange [*:*]
set cbtics autofreq scale 5.000000000000e-1,5.000000000000e-1
set border 15 front  lc rgb "black" lw 1
plot "-" binary endian=little record=26 format="%float64" using 1:2 with points pt 6 lc rgb "blue" t "S (Susceptible)", "-" binary endian=little record=26 format="%float64" using 1:2 with points pt 6 lc rgb "red" t "I (Infected)", "-" binary endian=little record=26 format="%float64" using 1:2 with points pt 6 lc rgb "green" t "R (Removed)"
             �Q@      �?�����lQ@       @(�r�Q@      @�`*�R@      @�R�t��R@      @�ܤN�(S@      @1�8G��S@      @f̄Ջ=T@       @�sp�T�T@      "@��5\3U@      $@x����U@      &@����S V@      (@�j,�WV@      *@14sk��V@      ,@��P��V@      .@�
<m�(W@      0@�a�c�^W@      1@҇�m
�W@      2@ �8|�W@      3@��lI��W@      4@Ə�*��W@      5@�|~
X@      6@��6X@      7@|sZMX@      8@�l��VaX@      9@5�� �rX@              >@      �?�����L8@       @�g��?�3@      @�k��]�/@      @>�׆u*@      @�)��0{%@      @R�"}��!@      @<)� N�@       @*:�#@      "@�>3W@      $@0N��;@      &@j�Og[W@      (@�-�M@      *@��3��@      ,@W�D3G~@      .@r{Z���@      0@B���!|�?      1@͛�����?      2@��q�(x�?      3@,�fޥ�?      4@HH�<�0�?      5@����\
�?      6@|�݌nO�?      7@�A�HF��?      8@���Ha�?      9@�}�k{�?                      �?������@       @L���Q#@      @�k$	R'@      @2o��
)@      @4�5��@)@      @&�I�(@      @2�S�)'@       @{W�B��%@      "@�{8��#@      $@(�����!@      &@��wI�' @      (@̉Gf6�@      *@)�=���@      ,@��M�@      .@w�,_�@      0@��V�1@      1@�k�5!@      2@�"�gd�@      3@��[��O	@      4@c���i@      5@�	�ځ�@      6@��v&q�@      7@��m@*�?      8@v�n����?      9@���J
��?
