#include "lnaimpedance.h"

const LNAImpedance::ctype LNAImpedance::impedanceArray[452] = {
	// 50 - 69 MHz.
	ctype(52.803,151.76), ctype(52.446,154.04), 
	ctype(52.19,157.03), ctype(51.811,160.69), 
	ctype(51.528,164.85), ctype(51.802,169.38), 
	ctype(52.157,174.17), ctype(52.522,179.32), 
	ctype(53.377,184.66), ctype(54.249,189.98), 
	ctype(55.512,195.65), ctype(56.949,201.94), 
	ctype(58.519,208.08), ctype(60.179,214.4), 
	ctype(62.755,220.82), ctype(65.071,227.42), 
	ctype(67.949,234.58), ctype(71.302,241.31), 
	ctype(74.795,248.26), ctype(78.823,255.1), 
	// 70 - 89 MHz.
	ctype(83.465,262.45), ctype(88.307,268.87), 
	ctype(93.864,275.2), ctype(99.639,281.3), 
	ctype(105.58,287.5), ctype(111.82,293.16), 
	ctype(118.57,298.85), ctype(125.37,303.59), 
	ctype(132.08,308.64), ctype(138.4,313.49), 
	ctype(144.44,318.07), ctype(150.64,323.03), 
	ctype(157.06,327.82), ctype(163.72,332.52), 
	ctype(170.03,337.71), ctype(176.67,343.27), 
	ctype(184.35,349.36), ctype(192.29,355.09), 
	ctype(200.02,360.97), ctype(208.58,366.41), 
	// 90 - 109 MHz.
	ctype(216.59,372.96), ctype(225.67,379.45), 
	ctype(235.53,386.17), ctype(245.42,391.84), 
	ctype(256.16,398.14), ctype(268.37,404.17), 
	ctype(281.23,410.06), ctype(295.19,416.49), 
	ctype(308.93,420.98), ctype(323.64,426.09), 
	ctype(339.04,431.74), ctype(354.61,436.51), 
	ctype(371.38,440.57), ctype(388.71,445.18), 
	ctype(406.16,447.12), ctype(426.15,450.65), 
	ctype(446.62,452.77), ctype(467.8,453.18), 
	ctype(489.21,450.95), ctype(510.8,449.09), 
	// 110 - 129 MHz.
	ctype(534.27,446.18), ctype(557.91,441.95), 
	ctype(582.01,436.09), ctype(605.65,427.46), 
	ctype(629.14,418.13), ctype(656.66,407.54), 
	ctype(681.79,394.13), ctype(705.83,378.75), 
	ctype(730.36,361.53), ctype(752.83,341.19), 
	ctype(775.31,320.56), ctype(798.29,296.29), 
	ctype(816.81,272.67), ctype(834.34,243.57), 
	ctype(851.51,217.39), ctype(865.02,186.63), 
	ctype(876.49,157.17), ctype(884.52,125.94), 
	ctype(891.00,96.191), ctype(894.14,61.696), 
	// 130 - 149 MHz.
	ctype(897.08,27.864), ctype(894.52,3.5946), 
	ctype(892.67,35.799), ctype(887.2,67.153), 
	ctype(881.6,99.435), ctype(872.62,129.87), 
	ctype(864.69,155.42), ctype(853.69,179.82), 
	ctype(840.15,205.61), ctype(825.47,227.86), 
	ctype(811.32,249.27), ctype(794.84,268.09), 
	ctype(778.56,283.31), ctype(760.98,301.41), 
	ctype(745.14,317.45), ctype(728.93,332.33), 
	ctype(712.62,344.89), ctype(695.47,357), 
	ctype(678.85,367.29), ctype(662.72,377.19), 
	// 150 - 169 MHz.
	ctype(645.2,384.83), ctype(628.67,392.63), 
	ctype(611.52,398.84), ctype(594.87,405.96), 
	ctype(577.85,411.64), ctype(561.2,416.75), 
	ctype(545.01,421.23), ctype(531.35,426.39), 
	ctype(516.58,428.47), ctype(501.43,431.9), 
	ctype(486.88,433.12), ctype(473.71,434.81), 
	ctype(461.38,435.33), ctype(448.13,435.89), 
	ctype(434.9,434.92), ctype(421.63,434.71), 
	ctype(409.75,433.7), ctype(397.58,432.03), 
	ctype(386.65,430.36), ctype(376.25,428.97), 
	// 170 - 189 MHz.
	ctype(366.39,427.51), ctype(356.49,426.02), 
	ctype(347.71,424.14), ctype(339.13,421.65), 
	ctype(330.76,420.74), ctype(322.54,417.59), 
	ctype(314.44,415.14), ctype(305.26,411.58), 
	ctype(297.84,408.67), ctype(290.45,405.71), 
	ctype(282.98,402.77), ctype(276.15,399.56), 
	ctype(268.95,397.06), ctype(262.05,393.77), 
	ctype(256.48,391.53), ctype(250.55,388.55), 
	ctype(244.67,385.83), ctype(239.24,383.2), 
	ctype(234.1,379.97), ctype(229.17,377.32), 
	// 190 - 209 MHz.
	ctype(224.07,374.82), ctype(219.45,372.16), 
	ctype(214.22,369.4), ctype(209.86,366.82), 
	ctype(205.51,363.94), ctype(201.02,361.7), 
	ctype(196.64,358.88), ctype(192.06,356.2), 
	ctype(187.99,353.37), ctype(184.01,350.71), 
	ctype(180.82,347.86), ctype(177.31,345.16), 
	ctype(173.86,342.35), ctype(170.04,339.97), 
	ctype(167.37,337.48), ctype(163.77,334.91), 
	ctype(160.71,332.13), ctype(157.16,329.43), 
	ctype(153.89,326.64), ctype(150.9,324.17), 
	// 210 - 229 MHz.
	ctype(148.42,321.22), ctype(145.32,318.78), 
	ctype(142.8,316.5), ctype(140.51,313.88), 
	ctype(138.07,311.21), ctype(135.65,308.95), 
	ctype(133.56,306.19), ctype(131.42,303.79), 
	ctype(129.42,301.32), ctype(127.57,298.73), 
	ctype(125.73,296.62), ctype(123.7,294.57), 
	ctype(121.82,292.36), ctype(119.93,290.17), 
	ctype(118.27,288.21), ctype(116.39,286.06), 
	ctype(114.48,283.88), ctype(112.76,281.88), 
	ctype(111.13,280.13), ctype(109.49,278.22), 
	// 230 - 249 MHz.
	ctype(107.74,276.41), ctype(105.99,274.7), 
	ctype(104.74,272.75), ctype(103.35,270.97), 
	ctype(101.76,269.11), ctype(100.21,267.31), 
	ctype(98.957,265.44), ctype(97.601,263.75), 
	ctype(96.383,261.71), ctype(94.926,260.1), 
	ctype(93.655,258.33), ctype(92.223,256.45), 
	ctype(91.028,254.56), ctype(89.743,252.72), 
	ctype(88.719,251.07), ctype(87.393,249.27), 
	ctype(86.092,247.57), ctype(84.946,245.82), 
	ctype(83.73,244.29), ctype(82.535,242.85), 
	// 250 - 269 MHz.
	ctype(81.653,241), ctype(80.563,239.12), 
	ctype(79.475,237.62), ctype(78.791,235.9), 
	ctype(77.636,234.28), ctype(76.866,232.99), 
	ctype(76.39,231.17), ctype(75.255,230.01), 
	ctype(74.317,228.71), ctype(73.493,227.35), 
	ctype(72.697,225.99), ctype(72.096,224.74), 
	ctype(71.661,223.36), ctype(70.576,222.31), 
	ctype(69.961,221.11), ctype(69.4,220.12), 
	ctype(68.796,218.83), ctype(67.94,217.5), 
	ctype(67.31,216.31), ctype(66.216,215.06), 
	// 270 - 289 MHz.
	ctype(65.584,213.54), ctype(64.877,212.28), 
	ctype(63.893,210.78), ctype(63.306,209.5), 
	ctype(62.675,208.11), ctype(61.585,207.05), 
	ctype(61.274,205.98), ctype(60.893,204.93), 
	ctype(60.162,203.75), ctype(59.491,202.39), 
	ctype(58.785,201.23), ctype(58.037,200.33), 
	ctype(57.763,199.1), ctype(57.082,197.79), 
	ctype(56.345,196.7), ctype(55.75,195.54), 
	ctype(55.49,194.62), ctype(54.995,193.59), 
	ctype(54.358,192.45), ctype(53.946,191.27), 
	// 290 - 309 MHz.
	ctype(53.429,190.25), ctype(52.941,189.14), 
	ctype(52.774,188.19), ctype(52.502,187.32), 
	ctype(52.019,186.36), ctype(51.827,185.47), 
	ctype(51.392,184.55), ctype(50.93,183.64), 
	ctype(50.507,182.67), ctype(50.048,181.69), 
	ctype(49.391,180.64), ctype(48.997,179.7), 
	ctype(48.551,178.69), ctype(48.098,177.85), 
	ctype(47.737,176.85), ctype(47.302,175.93), 
	ctype(46.857,174.99), ctype(46.454,174.01), 
	ctype(46.187,173.1), ctype(45.907,172.37), 
	// 310 - 329 MHz.
	ctype(45.59,171.47), ctype(45.215,170.71), 
	ctype(44.848,170.01), ctype(44.296,169.24), 
	ctype(43.968,168.56), ctype(43.492,167.7), 
	ctype(43.203,166.91), ctype(42.812,166.24), 
	ctype(42.453,165.32), ctype(42.256,164.43), 
	ctype(42.245,163.79), ctype(41.847,163.03), 
	ctype(41.537,162.33), ctype(41.211,161.44), 
	ctype(40.928,160.64), ctype(40.545,159.95), 
	ctype(40.022,159.18), ctype(39.594,158.33), 
	ctype(39.436,157.57), ctype(39.083,156.85), 
	// 330 - 349 MHz.
	ctype(38.585,155.97), ctype(38.308,155.22), 
	ctype(38.104,154.41), ctype(37.944,153.69), 
	ctype(37.686,152.85), ctype(37.312,152.05), 
	ctype(37.102,151.4), ctype(37.074,150.78), 
	ctype(36.776,149.94), ctype(36.544,149.26), 
	ctype(36.359,148.6), ctype(36.08,147.97), 
	ctype(35.99,147.3), ctype(35.746,146.63), 
	ctype(35.507,145.93), ctype(35.38,145.39), 
	ctype(35.17,144.62), ctype(34.934,143.98), 
	ctype(34.799,143.29), ctype(34.627,142.76), 
	// 350 - 369 MHz.
	ctype(34.442,142.15), ctype(34.097,141.54), 
	ctype(33.835,140.87), ctype(33.512,140.46), 
	ctype(33.297,139.72), ctype(32.991,139.29), 
	ctype(32.688,138.53), ctype(32.479,137.92), 
	ctype(32.298,137.37), ctype(32.054,136.76), 
	ctype(31.97,136.07), ctype(31.669,135.6), 
	ctype(31.519,134.87), ctype(31.281,134.45), 
	ctype(31.115,133.84), ctype(30.888,133.34), 
	ctype(30.746,132.8), ctype(30.581,132.2), 
	ctype(30.414,131.7), ctype(30.129,131.1), 
	// 370 - 389 MHz.
	ctype(30.007,130.48), ctype(29.826,129.88), 
	ctype(29.702,129.29), ctype(29.507,128.72), 
	ctype(29.274,128.22), ctype(29.115,127.6), 
	ctype(28.968,127.13), ctype(28.758,126.62), 
	ctype(28.543,126.12), ctype(28.386,125.6), 
	ctype(28.211,125.14), ctype(28.065,124.65), 
	ctype(27.857,124.23), ctype(27.758,123.82), 
	ctype(27.619,123.41), ctype(27.564,122.98), 
	ctype(27.532,122.53), ctype(27.377,122.06), 
	ctype(27.215,121.69), ctype(27.072,121.15), 
	// 390 - 409 MHz.
	ctype(26.882,120.59), ctype(26.655,120.07), 
	ctype(26.442,119.52), ctype(26.172,119.02), 
	ctype(25.979,118.44), ctype(25.805,117.82), 
	ctype(25.732,117.32), ctype(25.542,116.85), 
	ctype(25.435,116.32), ctype(25.251,115.89), 
	ctype(24.992,115.5), ctype(24.86,115.1), 
	ctype(24.668,114.68), ctype(24.544,114.24), 
	ctype(24.343,113.81), ctype(24.136,113.33), 
	ctype(23.995,112.84), ctype(23.9,112.27), 
	ctype(23.732,111.82), ctype(23.654,111.39), 
	// 410 - 429 MHz.
	ctype(23.413,110.92), ctype(23.39,110.41), 
	ctype(23.214,109.92), ctype(23.097,109.4), 
	ctype(22.958,108.96), ctype(22.836,108.51), 
	ctype(22.7,108.03), ctype(22.618,107.6), 
	ctype(22.428,107.26), ctype(22.419,106.93), 
	ctype(22.344,106.53), ctype(22.267,106.2), 
	ctype(22.137,105.8), ctype(22.05,105.3), 
	ctype(21.98,104.92), ctype(21.851,104.54), 
	ctype(21.762,104.11), ctype(21.648,103.81), 
	ctype(21.493,103.38), ctype(21.433,102.92), 
	// 430 - 449 MHz.
	ctype(21.34,102.64), ctype(21.154,102.32), 
	ctype(20.983,101.8), ctype(20.861,101.37), 
	ctype(20.606,101.01), ctype(20.523,100.61), 
	ctype(20.314,100.25), ctype(20.112,99.741), 
	ctype(19.995,99.266), ctype(19.955,98.997), 
	ctype(19.758,98.615), ctype(19.688,98.133), 
	ctype(19.582,97.774), ctype(19.581,97.42), 
	ctype(19.475,97.13), ctype(19.387,96.75), 
	ctype(19.285,96.317), ctype(19.257,95.976), 
	ctype(19.19,95.605), ctype(19.042,95.164), 
	// 450 - 469 MHz.
	ctype(18.918,94.76), ctype(18.813,94.341), 
	ctype(18.675,93.897), ctype(18.573,93.465), 
	ctype(18.45,93.078), ctype(18.416,92.723), 
	ctype(18.373,92.336), ctype(18.222,91.964), 
	ctype(18.151,91.644), ctype(18.142,91.32), 
	ctype(18.12,90.964), ctype(18.066,90.571), 
	ctype(17.921,90.204), ctype(17.831,89.837), 
	ctype(17.799,89.529), ctype(17.705,89.176), 
	ctype(17.555,88.883), ctype(17.417,88.52), 
	ctype(17.354,88.098), ctype(17.344,87.786), 
	// 470 - 489 MHz.
	ctype(17.25,87.435), ctype(17.144,87.005), 
	ctype(17.013,86.724), ctype(16.952,86.38), 
	ctype(16.872,86.119), ctype(16.786,85.864), 
	ctype(16.636,85.441), ctype(16.569,85.112), 
	ctype(16.519,84.808), ctype(16.461,84.443), 
	ctype(16.327,84.082), ctype(16.251,83.738), 
	ctype(16.066,83.371), ctype(15.943,83.089), 
	ctype(15.814,82.782), ctype(15.669,82.459), 
	ctype(15.522,82.064), ctype(15.458,81.733), 
	ctype(15.293,81.351), ctype(15.243,81.042), 
	// 490 - 500 MHz.
	ctype(15.146,80.671), ctype(15.05,80.296), 
	ctype(14.973,79.911), ctype(14.955,79.614), 
	ctype(14.876,79.266), ctype(14.839,78.892), 
	ctype(14.723,78.569), ctype(14.683,78.282), 
	ctype(14.643,78.058), ctype(14.614,77.894), 
	ctype(14.594,77.713),
	// Extra element to allow 'interpolation' up to 500 MHz
	ctype(14.594,77.713)
};

