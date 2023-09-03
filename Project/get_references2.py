import pandas as pd
import subprocess
from itertools import product
import os
import random
import time
import re

def extend_ambiguous_dna(seq):
	degenerated_table = {
		'A': ['A'],
		'C': ['C'],
		'G': ['G'],
		'T': ['T'],
		'W': ['A', 'T'],
		'S': ['C', 'G'],
		'M': ['A', 'C'],
		'K': ['G', 'T'],
		'R': ['A', 'G'],
		'Y': ['C', 'T'],
		'B': ['C', 'G', 'T'],
		'D': ['A', 'G', 'T'],
		'H': ['A', 'C', 'T'],
		'V': ['A', 'C', 'G'],
		'N': ['A', 'C', 'G', 'T'],
	}

	new_seq = ''
	for base in seq:
		if base == 'I':
			new_seq += 'N'
		else:
			new_seq += base
	seq = new_seq

	r = []
	# Aplica produto cartesiano nos conjuntos (conjunto = possíveis bases para cada letra)
	for i in product(*[degenerated_table[j] for j in seq]):
		r.append("".join(i))

	return r


def create_fasta(individual_seqs, output_file):
	with open(output_file, "w") as fasta_file:
		for primer_code, seq in individual_seqs:
			fasta_record = f">{primer_code}\n{seq}\n"
			fasta_file.write(fasta_record)


def run_blast(fasta_file, output_file, max_target_seqs='10', num_threads='6', database='/media/bioinfo/6tb_hdd/04_Blast_Databases/BLAST_DB_nt/nt'):
	subprocess.run(
		f'export BLASTDB=/media/bioinfo/6tb_hdd/04_Blast_Databases/taxdb && \
		time blastn \
		-db {database} \
		-query {fasta_file} \
		-outfmt "6 qseqid pident gaps qcovs mismatch evalue score bitscore sgi sacc slen stitle sscinames" \
		-out {output_file} \
		-num_threads {num_threads} \
		-max_target_seqs {max_target_seqs} \
		-task blastn-short',
		shell=True,
		executable='/bin/bash'
	)


def split_and_duplicate(row):
    if '@' in row:
        valores_divididos = row.split('@')
        primer = ''.join(valores_divididos[0:-2])
        pos = valores_divididos[-1]
        return pd.Series([valores_divididos[0], pos])
    else:
        return pd.Series([row, row])


def download_references(df_references, primer_code, max_len, count, selected_organisms, path2save_references):
	df_primer_references = df_references.loc[df_references['Primer'] == primer_code].reset_index(drop=True)

	if max_len is not None:
		df_primer_references = df_primer_references.loc[df_primer_references['slen'] <= max_len].reset_index(drop=True)
		df_primer_references = df_primer_references.loc[df_primer_references['pident'] >= 95].reset_index(drop=True)
		df_primer_references = df_primer_references.loc[df_primer_references['qcovs'] >= 95].reset_index(drop=True)

	df_primer_references = df_primer_references.drop_duplicates(subset=['sscinames'])
	df_primer_references = df_primer_references.dropna(subset=['sscinames'])

	for organism, sgi, db in zip(df_primer_references['sscinames'], df_primer_references['sgi'], df_primer_references['level_0']):
		organism_renamed = organism.replace(' ', '_')
		try:
			if organism not in selected_organisms:

				if db == 'nt':
					subprocess.run(
						f'cd {path2save_references}{primer_code}/ &&\
						blastdbcmd -db /media/bioinfo/6tb_hdd/04_Blast_Databases/BLAST_DB_nt/nt -entry {sgi} > {organism_renamed}.fasta',
						shell=True,
						executable='/bin/bash'
					)
					count += 1
					selected_organisms.append(organism)
				elif db == '16S':
					subprocess.run(
						f'cd {path2save_references}{primer_code}/ &&\
						blastdbcmd -db /media/bioinfo/6tb_hdd/04_Blast_Databases/16sDB/16S_ribosomal_RNA -entry {sgi} > {organism_renamed}.fasta',
						shell=True,
						executable='/bin/bash'
					)
					selected_organisms.append(organism)
					count += 1

		except Exception as e:
			print(e)

		if count == 10:
			return count, selected_organisms

	return count, selected_organisms

if __name__ == '__main__':
	path2save_references = '/media/bioinfo/6tb_hdd/03_ELLEN/02_data/primers_references/'
	primers = pd.read_csv('primers.csv', sep=';')
	
	filtered_primers = [3159, 1184, 765, 2404, 2406, 2405, 2407, 2414, 2415, 2418, 1088, 301, 3420, 192, 193, 189, 191, 757, 758, 1443, 960, 88, 3268, 1823, 1824, 1013, 1016, 4, 3703, 690, 3275, 3648, 3679, 3649, 3680, 3787, 3755, 3756, 3792, 3770, 3748, 3743, 3757, 3769, 3788, 3750, 3759, 79, 232, 1787, 1789, 1786, 1788, 1790, 216, 974, 217, 218, 622, 1631, 1632, 1633, 1634, 1635, 1636, 1637, 1638, 1639, 1640, 2025, 2028, 863, 864, 3235, 2986, 22, 23, 24, 25, 2277, 2526, 2524, 2527, 2525, 2475, 772, 773, 3055, 1203, 1427, 1428, 3697, 2814, 3404, 876, 877, 2318, 3043, 878, 879, 880, 2291, 875, 874, 1984, 2872, 2319, 3690, 1104, 1357, 1370, 82, 3510, 2607, 412, 413, 414, 1198, 1199, 3088, 3722, 1658, 2795, 220, 221, 2960, 534, 1910, 3185, 563, 3631, 2752, 212, 213, 277, 811, 276, 2724, 2518, 2528, 381, 385, 387, 389, 380, 384, 388, 222, 70, 71, 3151, 1419, 1521, 683, 842, 841, 820, 225, 1691, 3172, 390, 1348, 2977, 1349, 1528, 2979, 2981, 3080, 2644, 3771, 367, 1750, 3785, 2995, 366, 3724, 1692, 391, 1530, 1531, 3167, 3166, 2688, 1586, 3413, 3414, 3153, 2144, 1398, 1686, 2084, 2138, 2139, 243, 1812, 608, 1985, 3090, 319, 318, 3376, 317, 242, 2462, 2461, 3716, 3717, 298, 299, 360, 533, 2520, 2876, 2059, 2622, 3145, 200, 199, 1344, 1345, 175, 3733, 3067, 3066, 1580, 3721, 615, 2250, 2251, 227, 265, 1576, 1577, 2182, 748, 60, 747, 749, 2422, 2423, 1579, 3391, 3394, 3393, 3392, 2975, 80, 2012, 2011, 2013, 3708, 61, 196, 197, 198, 3440, 3442, 124, 1717, 581, 3033, 73, 72, 3602, 254, 255, 20, 21, 3695, 1020, 3328, 62, 324, 166, 170, 172, 187, 983, 174, 607, 1794, 2897, 1974, 3082, 2493, 336, 647, 1102, 554, 846, 3214, 2900, 2246, 3122, 59, 840, 1410, 2483, 3371, 912, 2244, 3232, 1078, 3430, 2430, 3684, 1784, 3599, 3277, 2397, 2522, 1317, 1318, 3397, 3398, 3399, 2410, 2411, 665, 1472, 416, 2581, 2974, 763, 3060, 409, 1101, 2899, 1258, 283, 3725, 2247, 58, 839, 893, 1377, 427, 1334, 3786, 2875, 2, 1666, 1877, 207, 866, 304, 307, 2074, 1681, 1234, 3692, 3693, 3142, 3395, 3396, 115, 1800, 1695, 2193, 2412, 2413, 2416, 2417, 2420, 2419, 2686, 2687, 524, 2599, 2600, 2598, 2188, 3329, 995, 55, 845, 751, 2517, 3630, 401, 745, 746, 3734, 2204, 2203, 1573, 1116, 2038, 2039, 1092, 1093, 392, 393, 374, 2989, 2990, 2991, 57, 1803, 1832, 291, 195, 194, 257, 2608, 2609, 2610, 1095, 1949, 1956, 1957, 1958, 1959, 1960, 1961, 1950, 1951, 1952, 1953, 1955, 1965, 1966, 2920, 281, 2916, 2921, 2923, 2917, 2918, 2927, 2922, 424, 869, 1725, 1356, 1698, 1699, 1700, 583, 1265, 623, 844, 1021, 1022, 2007, 3709, 2069, 2070, 1154, 2067, 2068, 1380, 967, 1468, 3441, 1024, 1023, 1826, 87, 16, 2721, 2579, 2557, 1025, 1332, 582, 584, 1270, 1371, 1827, 3661, 3662, 2481, 2482, 3084, 3085, 2972, 2973, 1391, 827, 1651, 1873, 1872, 1321, 2425, 2611, 1322, 2400, 2401, 2760, 766, 2558, 1876, 1151, 1152, 644, 643, 2746, 896, 1687, 1202, 549, 2268, 1816, 1144, 1145, 1147, 320, 321, 2195, 3735, 3028, 1113, 1114, 5, 8, 18, 7, 9, 3409, 3408, 3407, 3410, 1175, 1176, 1393, 1394, 3182, 3379, 2130, 2647, 2129, 1837, 629, 1836, 632, 2369, 2370, 3363, 3618, 3617, 1962, 2934, 2136, 2769, 3236, 3194, 2770, 3237, 3195, 2850, 2398, 1065, 2684, 1066, 578, 2399, 1064, 631, 2650, 2292, 2293, 2508, 2821, 2033, 3590, 3300, 2032, 3388, 2031, 3252, 3587, 1408, 3588, 3223, 2820, 547, 2176, 2712, 3311, 3351, 3383, 2433, 600, 1730, 3729, 3730, 2815, 3384, 3188, 661, 3056, 2432, 1298, 1299, 1885, 1493, 1494, 2932, 2933, 2045, 2046, 2049, 2050, 3313, 1401, 3302, 1430, 2377, 1923, 2043, 2044, 1359, 1360, 1302, 1308, 1306, 1303, 1307, 2714, 1304, 1300, 717, 716, 606, 2445, 2235, 2895, 2391, 2125, 2126, 2447, 3715, 2127, 1741, 2448, 2393, 2446, 3714, 3119, 1742, 1711, 3200, 1759, 2541, 2540, 2894, 790, 1758, 2395, 1470, 2304, 3738, 3668, 2305, 2209, 3712, 2271, 3765, 3332, 3713, 3490, 2683, 2574, 2272, 3766, 1469, 2164, 3149, 1749, 2256, 1672, 2257, 3706, 3707, 519, 436, 447, 479, 514, 448, 462, 3234, 292, 2839, 3623, 3624, 2567, 2569, 1527, 3190, 3191, 3323, 3324, 2234, 2165, 1499, 1500, 2253, 1675, 1501, 3233, 3745, 3746, 2206, 2207, 792, 2208, 3326, 3325, 1389, 1390, 2308, 667, 668, 664, 1294, 871, 3322, 3330, 2362, 2380, 2381, 2506, 2236, 2505, 2237, 3793, 2913, 571, 572, 420, 1630, 1295, 3698, 2503, 2325, 1629, 3331, 1678, 1679, 3189, 2504, 2326, 2623, 2617, 2450, 2124, 2350, 2351, 2628, 953, 956, 3468, 1539, 363, 1821, 364, 1542, 1376, 3499, 1943, 3034, 1405, 954, 755, 2361, 2258, 3767, 1434, 1436, 3699, 663, 3768, 756, 870, 1763, 1765, 1767, 1762, 3320, 1768, 1769, 1770, 1772, 1773, 3660, 3659, 1073, 2620, 558, 3063, 786, 788, 3030, 3031, 2757, 1932, 2424, 1781, 1707, 1207, 2054, 613, 2055, 1087, 1075, 1076, 3156, 3157, 3154, 2621, 3155, 695, 696, 1090, 1091, 2178, 1808, 1809, 2230, 1206, 1706, 1766, 1764, 3039, 3367, 1190, 3041, 3507, 910, 977, 861, 552, 2167, 2168, 2835, 1480, 3270, 3495, 3040, 2373, 589, 3271, 3615, 2320, 968, 969, 587, 2640, 2992, 1774, 3657, 3658, 722, 739, 723, 725, 729, 735, 740, 727, 734, 1228, 1215, 1217, 3632, 1455, 348, 1908, 2468, 1196, 2490, 609, 345, 430, 672, 1197, 344, 1454, 1907, 1685, 1079, 812, 3202, 1585, 1259, 2467, 2345, 2222, 2221, 3651, 3652, 1226, 738, 728, 732, 720, 724, 726, 730, 736, 741, 733, 721, 1222, 3650, 355, 2914, 2509, 3021, 356, 3046, 1924, 685, 2366, 83, 2212, 3783, 3784, 2211, 753, 536, 1888, 1235, 857, 858, 2347, 2348, 2491, 2492, 798, 84, 1584, 899, 1342, 2854, 997, 1335, 3365, 85, 754, 2853, 1343, 2855, 2488, 311, 309, 313, 527, 3165, 3366, 310, 3181, 2856, 2857, 1320, 1886, 3672, 684, 774, 771, 3068, 3069, 3671, 1891, 397, 1904, 395, 396, 3576, 3578, 3580, 3577, 3579, 3581, 279, 1921, 278, 2744, 2562, 898, 2785, 638, 637, 1388, 2587, 796, 1245, 1246, 1248, 1249, 2560, 2281, 2282, 814, 550, 406, 3494, 1387, 405, 404, 628, 354, 538, 400, 2371, 2924, 340, 979, 980, 338, 339, 3171, 1869, 3796, 1482, 1483, 2477, 2478, 1565, 2057, 907, 2937, 2938, 2939, 3083, 2363, 3253, 1063, 708, 835, 836, 1059, 707, 2852, 2851, 2711, 2364, 1548, 2612, 2629, 2630, 828, 2789, 709, 3640, 3432, 966, 2985, 830, 2586, 1866, 1867, 1865, 3434, 3436, 3438, 3435, 3437, 3439, 837, 965, 2984, 1733, 2585, 1099, 1100, 692, 2135, 3433, 1563, 3622, 1655, 1656, 833, 829, 3641, 710, 349, 3663, 2830, 1552, 1553, 2149, 1186, 1187, 2210, 3222, 275, 742, 743, 274, 1543, 2546, 2829, 1352, 1351, 625, 2864, 1268, 1269, 3373, 541, 542, 2865, 544, 656, 651, 657, 3251, 3634, 3250, 1465, 2501, 2140, 1466, 601, 602, 603, 2196, 2197, 604, 2502, 2141, 3656, 2386, 2387, 2388, 1001, 1002, 652, 3192, 3193, 660, 659, 658, 654, 3592, 3593, 3751, 3752, 947, 1097, 1098, 1291, 2822, 1293, 2532, 1729, 1726, 689, 2147, 2148, 3354, 3353, 3719, 3484, 2232, 2233, 1296, 2902, 2903, 2870, 2901, 2592, 2594, 2871, 1727, 1728, 1327, 1329, 1328, 1069, 1285, 1286, 3053, 1287, 3261, 1134, 1870, 1135, 1806, 889, 825, 824, 826, 1007, 823, 1309, 1310, 890, 1549, 1138, 1137, 1139, 1136, 1484, 1485, 1311, 1312, 3012, 3013, 2512, 2514, 944, 2511, 2513, 945, 3619, 3620, 2862, 2863, 885, 888, 886, 887, 1325, 1326, 701, 702, 700, 981, 3421, 3422, 697, 698, 567, 624, 564, 565, 3799, 3798, 2515, 1694, 3029, 1708, 2649, 598, 1945, 2798, 1796, 1798, 1797, 1649, 1650, 3146, 3603, 1785, 3118, 1557, 3036, 3037, 2539, 2538, 2040, 2041, 1058, 1833, 1813, 1055, 2266, 2632, 1947, 1948, 2619, 1880, 1881, 2131, 2132, 1533, 3369, 1323, 2352, 2353, 1324, 1857, 1858, 3196, 648, 649, 3197, 1284, 3212, 3772, 868, 3011, 867, 3010, 3238, 2910, 2911, 999, 2838, 1156, 1155, 2837, 676, 675, 2331, 706, 705, 872, 873, 3591, 3305, 1036, 1037, 2486, 925, 931, 3049, 3050, 1056, 1057, 2479, 2480, 3242, 3244, 3246, 3243, 3266, 3267, 3509, 1251, 2466, 3508, 2836, 1252, 1859, 1385, 3710, 2970, 2971, 1281, 1282, 2826, 2199, 2827, 2828, 1735, 1736, 2614, 2615, 1278, 1280, 2223, 1439, 1442, 3676, 2928, 2929, 3024, 3025, 3092, 1599, 3654, 1532, 1597, 1594, 3003, 1421, 1420, 1739, 1747, 1536, 2523, 2324, 3163, 3162, 1663, 3318, 1554, 3319, 3506, 2228, 2229, 3262, 2214, 3309, 3160, 3161, 3497, 3498, 1710, 3340, 3341, 1917, 1918, 1920, 2536, 2537, 2020, 2021, 2022, 2023, 2375, 3179, 2500, 3203, 2255, 2860, 2861, 2858, 2793, 2794, 2455, 3669, 2530, 2531, 3675, 3673, 2763, 3306, 3415, 3416, 3307, 2867, 2866, 3124, 3125, 3126, 3127, 3128, 3129, 3130, 3131, 3132, 3133, 3134, 3135, 3136, 3137, 3138, 3139, 3140, 3447, 3272, 3273, 3411, 3412, 3472, 3471, 3469, 3470, 3762, 3763, 3790, 3760, 3761, 3773, 3774, 3775, 3776, 3777, 3778]


	redo = [2291, 2960, 390, 3080, 391, 3413, 3414, 608, 1985, 3090, 2250, 3708, 1717, 554, 1695, 2609, 2007, 3709, 2425, 320, 321, 2129, 661, 1298, 1299, 3302, 1359, 1360, 2125, 2126, 2127, 1711, 790, 3490, 2567, 2569, 2380, 420, 2503, 2504, 2623, 2617, 2467, 536, 84, 1891, 898, 2371, 836, 2501, 2196, 2197, 2502, 3656, 2592, 1327, 1329, 1328, 1285, 1286, 1287, 1134, 1135, 1309, 1310, 1484, 1485, 1311, 1312, 2515, 1385, 2614, 2615, 3092, 3506, 1710, 2500, 3203, 3675, 2763, 3447]

	redo = [3709, 2007]

	
	seqs = []
	for primer_id, primer_code, primer_seq in zip(primers['id'], primers['code'], primers['nuc']):
		if primer_id in filtered_primers and primer_id in redo:
			primer_seq_extended = extend_ambiguous_dna(primer_seq)
			primer_code = primer_code.replace(' ', '_')
			primer_code = primer_code.replace('(', '-')
			primer_code = primer_code.replace(')', '-')
			primer_code = primer_code.replace('/', '-')
			full_pattern = re.compile('[^a-zA-Z0-9_-]')
			primer_code = str(re.sub(full_pattern, '_', primer_code))
			
			try:
				if len(primer_seq_extended) == 1:
					seqs.append((f'{primer_code}@{primer_code}', primer_seq_extended[0]))
				else:
					if len(primer_seq_extended) > 50:
						primer_seq_extended = random.choices(primer_seq_extended, k=50)  # talvez trocar pra 10% ou 25% da quantidade total inves de deixar fixo 50

					for pos, extended_primer in enumerate(primer_seq_extended):
						seqs.append((f'{primer_code}@{pos}', extended_primer))

			except:
				print(primer_code)

	create_fasta(seqs, f'get_references2/primers_seqs_redo2.fa')

	print('Rodando BLAST NT')
	run_blast(f'./get_references2/primers_seqs_redo2.fa', f'./get_references2/primers_seqs_redo2_nt.tsv', max_target_seqs='50', num_threads='12')
	print('Rodando BLAST 16S')
	run_blast(f'./get_references2/primers_seqs_redo2.fa', f'./get_references2/primers_seqs_redo2_16S.tsv', max_target_seqs='50', num_threads='12', database='/media/bioinfo/6tb_hdd/04_Blast_Databases/16sDB/16S_ribosomal_RNA')

	start = time.process_time()

	df_references_nt = pd.read_csv(f'get_references2/primers_seqs_redo2_nt.tsv', sep='\t', names=['qseqid', 'pident', 'gaps', 'qcovs', 'mismatch', 'evalue', 'score', 'bitscore', 'sgi', 'sacc', 'slen', 'stitle', 'sscinames'])
	df_references_16S = pd.read_csv(f'get_references2/primers_seqs_redo2_16S.tsv', sep='\t', names=['qseqid', 'pident', 'gaps', 'qcovs', 'mismatch', 'evalue', 'score', 'bitscore', 'sgi', 'sacc', 'slen', 'stitle', 'sscinames'])

	df_references = pd.concat([df_references_nt, df_references_16S], keys=['nt', '16S']).reset_index()
	df_references['sscinames'] = df_references['sscinames'].str.replace(' ', '_')
	df_references['sscinames'] = df_references['sscinames'].str.replace('(', '-')
	df_references['sscinames'] = df_references['sscinames'].str.replace(')', '-')

	full_pattern = re.compile('[^a-zA-Z0-9_-]')
	df_references['sscinames'] = df_references['sscinames'].str.replace(full_pattern, '_', regex=True)


	df_references[['Primer', 'SeqID']] = df_references['qseqid'].apply(split_and_duplicate)

	df_references = df_references.sample(frac = 1).reset_index(drop=True)

	df_references = df_references.sort_values(by=['qcovs', 'pident', 'slen'], ascending=[False, False, False]).reset_index(drop=True)

	primer_codes = df_references['Primer'].unique()

	# for primer_id, primer_code, primer_seq in zip(primers['id'], primers['code'], primers['nuc']):
	# 	primer_code = str(primer_code).replace('/', '_-_')
	# 	if primer_id in filtered_primers and primer_code not in primer_codes:
	# 		print(primer_code)

	# primers_blast = [str(primer_code) for primer_code in primer_codes]
	# primers_original = [str(primer_code) for primer_code in primers['code']]
	# for primer_code in primers_blast:
	# 	if primer_code not in primers_original:
	# 		print('diferente', primer_code)

	# print(primers_blast)
	# print(primers_original)

	redo_primer_codes = []
	for primer_id, primer_code in zip(primers['id'], primers['code']):
		if primer_id in redo:
			primer_code = primer_code.replace(' ', '_')
			primer_code = primer_code.replace('(', '-')
			primer_code = primer_code.replace(')', '-')
			primer_code = primer_code.replace('/', '-')
			
			full_pattern = re.compile('[^a-zA-Z0-9_-]')
			primer_code = str(re.sub(full_pattern, '_', primer_code))
			redo_primer_codes.append(primer_code)

	for primer_code in primer_codes:
		if primer_code in redo_primer_codes:
			print(primer_code)

			if not os.path.exists(f"{path2save_references}{primer_code}"):
				os.makedirs(f"{path2save_references}{primer_code}")
				
				count = 0
				max_len = (30000, 40000, 50000, 100000, 200000)
				len_pos = 0
				selected_organisms = []
				while count < 10:
					if len_pos == 5:
						count, selected_organisms = download_references(df_references, primer_code, None, count, selected_organisms, path2save_references)
					elif len_pos == 6:
						break
					else:
						count, selected_organisms = download_references(df_references, primer_code, max_len[len_pos], count, selected_organisms, path2save_references)
					len_pos += 1

			# df_primer_references = df_references.loc[df_references['Primer'] == primer_code].reset_index(drop=True)

			# df_primer_references = df_primer_references.loc[df_primer_references['slen'] < 50000].reset_index(drop=True)

			# df_primer_references = df_primer_references.drop_duplicates(subset=['sscinames'])

			# count = 0
			# for organism, sgi in zip(df_primer_references['sscinames'], df_primer_references['sgi']):
			# 	organism_renamed = organism.replace(' ', '_')
			# 	try:
			# 		subprocess.run(
			# 			f'cd {path2save_references}{primer_code}/ &&\
			# 			blastdbcmd -db /media/bioinfo/6tb_hdd/04_Blast_Databases/BLAST_DB_nt/nt -entry {sgi} > {organism_renamed}.fasta',
			# 			shell=True,
			# 			executable='/bin/bash'
			# 		)
			# 		count += 1
			# 	except Exception as e:
			# 		print(e)

			# 	if count == 10:
			# 		break

	end = time.process_time()
	total = end - start

	print('Tempo para download de sequencias: ', total)
