import subprocess
import pandas as pd
from statistics import mode


def run_blast(max_target_seqs='12000', num_threads='10', database='/media/bioinfo/6tb_hdd/04_Blast_Databases/PrimersReferences/primers_references'):
	subprocess.run(
		f'export BLASTDB=/media/bioinfo/6tb_hdd/04_Blast_Databases/taxdb && \
		time blastn \
		-db /media/bioinfo/6tb_hdd/04_Blast_Databases/PrimersReferences/primers_references \
		-query primers.fasta \
		-outfmt "6 qseqid pident gaps qcovs mismatch evalue score bitscore sgi sacc slen stitle sscinames" \
		-out /media/bioinfo/6tb_hdd/03_ELLEN/02_data/primers_references_blast_alignment/blast_output.tsv \
		-num_threads 10 \
		-max_target_seqs 12000 \
		-task blastn-short',
		shell=True,
		executable='/bin/bash'
	)


def split_and_duplicate_primer(row, sep):
    if sep in row:
        splitted_values = row.split(sep)
        primer = ''.join(splitted_values[0:-2])
        pos = splitted_values[-1]
        return pd.Series([splitted_values[0], pos])
    else:
        return pd.Series([row, row])


def split_orientation(row, sep):
    if sep in row:
        splitted_values = row.split(sep)
        return splitted_values[-1]


def get_primers_alignment_data(primers, df):
	primers_alignment_data = {}
	for primer_id in primers:

		primer_df = df.loc[df['Primer'] == primer_id].reset_index(drop=True)

		primers_alignment_data[primer_id] = primer_df

	return primers_alignment_data


def clean_primers_alignment_data(primers_alignment_data, primers):
	for item in primers_alignment_data.items():
		primer_id = item[0]
		primer_df = item[1]
		primer_seq = primers.loc[primers['id'] == primer_id]['nuc'].iloc[0]
		
		# alinhamentos do tamanho do primer até o tamanho do primer - 2
		primer_df = primer_df.loc[primer_df['length'].isin((len(primer_seq), len(primer_seq) - 1, len(primer_seq) - 2))]

		# até 2 mismatches
		primer_df = primer_df.loc[primer_df['mismatch'] <= 2]

		primers_alignment_data[primer_id] = primer_df

	return primers_alignment_data


def get_common_organisms(alignments_0_df, alignments_1_df):		
	organisms_0 = set(alignments_0_df['sacc'])
	
	organisms_1 = set(alignments_1_df['sacc'])

	common_organisms = organisms_0.intersection(organisms_1)
	return common_organisms


def check_product_size(alignments_0_df, alignments_1_df, organism):
	organism_0 = alignments_0_df.loc[alignments_0_df['sacc'] == organism].reset_index(drop=True)
	pos_0 = organism_0['sstart'][0]

	organism_1 = alignments_1_df.loc[alignments_1_df['sacc'] == organism].reset_index(drop=True)
	pos_1 = organism_1['send'][0]

	if abs(pos_0 - pos_1) >= 200 and abs(pos_0 - pos_1) <= 10000:
		return True
	else:
		return False


def get_pair_orientation(alignments_0_df, alignments_1_df, pair, organism):
	organism_0 = alignments_0_df.loc[alignments_0_df['sacc'] == organism].reset_index(drop=True)
	pos_0 = organism_0['sstart'][0]
	orientation_0 = organism_0['orientation'][0]

	organism_1 = alignments_1_df.loc[alignments_1_df['sacc'] == organism].reset_index(drop=True)
	pos_1 = organism_1['sstart'][0]
	orientation_1 = organism_1['orientation'][0]

	orientation = None
	if pos_0 < pos_1 and orientation_0 != '-1' and orientation_1 != '-1' :
		orientation = (pair[0], pair[1])
	elif pos_0 > pos_1 and orientation_0 != '-1' and orientation_1 == '-1':
		orientation = (pair[0], pair[1])
	elif pos_0 < pos_1 and orientation_0 != '-1' and orientation_1 == '-1':
		orientation = (pair[1], pair[0])
	elif pos_0 > pos_1 and orientation_0 != '-1' and orientation_1 != '-1':
		orientation = (pair[1], pair[0])
	elif pos_0 < pos_1 and orientation_0 == '-1' and orientation_1 != '-1':
		orientation = (pair[1], pair[0])
	elif pos_0 > pos_1 and orientation_0 == '-1' and orientation_1 != '-1':
		orientation = (pair[0], pair[1])
	elif pos_0 < pos_1 and orientation_0 == '-1' and orientation_1 == '-1':
		orientation = (pair[0], pair[1])
	elif pos_0 > pos_1 and orientation_0 == '-1' and orientation_1 == '-1':
		orientation = (pair[1], pair[0])
	else:
		print(pos_0, pos_1, orientation_0, orientation_1)

	return orientation


def get_primers_pairs(primers_alignment_data):
	primer_pairs = [(primer_a, primer_b) for idx, primer_a in enumerate(primers_alignment_data.keys()) for primer_b in list(primers_alignment_data.keys())[idx + 1:]]

	print(f'Total de pares a serem processados: {len(primer_pairs)}')
	print('\n\n\n')
	accepted_pairs = []
	for pos, pair in enumerate(primer_pairs):
		common_organisms = get_common_organisms(primers_alignment_data[pair[0]], primers_alignment_data[pair[1]])

		orientations = []
		if len(common_organisms) >= 20:  # pelo menos metade das referencias são comun entre eles 
			for organism in common_organisms:
				if check_product_size(primers_alignment_data[pair[0]], primers_alignment_data[pair[1]], organism):
					pair_orientation = get_pair_orientation(primers_alignment_data[pair[0]], primers_alignment_data[pair[1]], pair, organism)
					if pair_orientation is not None:
						orientations.append(pair_orientation)

		if len(orientations) >= 20:
			print('\tPar aceito, orientação:', mode(orientations))
			accepted_pairs.append(mode(orientations))

		print(f'Par atual: {pos + 1}')

	return accepted_pairs


def write_pairs_file(accepted_pairs, primers_meta):
	pairs_dict = {'NAME': [], 'LEFT_PRIMER': [], 'RIGHT_PRIMER': []}
	for pair in accepted_pairs:
		name = f'{pair[0]}_{pair[1]}'
		LEFT_PRIMER = primers_meta.loc[primers_meta['id'] == pair[0]]['nuc'].iloc[0]  # ver que tem primers que eu renomeei no outro programa.... isso na verdade tem que ser uma limpeza inclusa no programa de extração dos primers do bold
		RIGHT_PRIMER = primers_meta.loc[primers_meta['id'] == pair[1]]['nuc'].iloc[0]

		pairs_dict['NAME'].append(name)
		pairs_dict['LEFT_PRIMER'].append(LEFT_PRIMER)
		pairs_dict['RIGHT_PRIMER'].append(RIGHT_PRIMER)

	pd.DataFrame(pairs_dict).to_csv('primers_pairs.tsv', sep='\t', index=False)


if __name__ == '__main__':
	filtered_primers = ['1136', '2044', '2318', '2387', '1687', '2851', '931', '2041', '1759', '2325', '88', '213', '3392', '2012', '3238', '2623', '427', '3031', '2617', '536', '1063', '725', '1215', '3785', '1832', '1401', '2149', '1557', '1352', '292', '1134', '2244', '2867', '1966', '820', '1888', '2266', '2045', '1750', '1781', '2910', '1298', '2407', '3729', '1300', '2445', '2560', '3172', '3275', '344', '58', '1349', '1579', '1837', '3185', '2414', '2821', '1308', '2585', '734', '79', '582', '1147', '1918', '2292', '1335', '3470', '72', '3188', '751', '2135', '2425', '702', '2611', '705', '2432', '741', '1679', '2826', '844', '1099', '2022', '648', '1749', '3699', '380', '385', '2038', '1097', '564', '3353', '1199', '2970', '1059', '1145', '2607', '2579', '3326', '3784', '3698', '638', '2540', '2933', '3088', '1468', '960', '628', '3167', '728', '3634', '3591', '2793', '1076', '3237', '3331', '309', '1634', '3712', '3028', '2461', '3721', '774', '1586', '2752', '2272', '1442', '3680', '945', '174', '3369', '2934', '717', '199', '1294', '3766', '3397', '2164', '864', '3160', '1543', '1228', '1747', '710', '1650', '1421', '3267', '1796', '2989', '55', '1554', '3300', '1536', '1885', '1345', '2380', '1786', '624', '2974', '2488', '514', '3734', '2481', '2557', '726', '2850', '1952', '899', '3234', '1376', '974', '2399', '20', '61', '1763', '124', '254', '629', '338', '1102', '2416', '1480', '868', '2199', '2126', '692', '2876', '2418', '3777', '647', '1184', '1469', '2503', '317', '1859', '2619', '2020', '1597', '3272', '1917', '1553', '414', '2230', '3735', '1762', '1098', '2124', '748', '2760', '3371', '1880', '2514', '3157', '3251', '997', '2050', '1638', '1427', '907', '2508', '2364', '3669', '1303', '1371', '1580', '667', '1631', '2820', '2031', '826', '1186', '283', '631', '1672', '3253', '1741', '2837', '745', '2517', '2450', '863', '3602', '544', '2405', '2712', '276', '3153', '194', '833', '1965', '1630', '3365', '572', '2518', '2351', '291', '3085', '690', '3197', '1959', '2025', '2195', '652', '393', '3684', '2412', '227', '875', '339', '3182', '1877', '722', '735', '3509', '3759', '3792', '2415', '2049', '217', '549', '874', '966', '404', '2852', '3714', '1824', '660', '643', '707', '672', '1736', '2899', '24', '3396', '2512', '547', '2629', '3243', '3703', '2433', '552', '3011', '603', '200', '416', '2829', '348', '1370', '2683', '242', '3577', '1787', '2211', '601', '1694', '1920', '1958', '3084', '395', '1876', '2423', '1710', '1501', '3151', '3203', '2381', '541', '1960', '3068', '2795', '3066', '2448', '311', '1322', '1548', '1025', '3769', '1858', '3439', '2822', '1770', '3395', '1728', '1955', '721', '2023', '1773', '1867', '2644', '2483', '281', '766', '1278', '1768', '3082', '2567', '3195', '396', '218', '83', '2830', '2902', '1296', '3743', '3366', '2866', '397', '2856', '2586', '1351', '3495', '2541', '2530', '3162', '2688', '983', '2486', '615', '366', '2501', '3469', '2527', '1635', '1789', '885', '1886', '2981', '772', '2221', '1332', '2509', '2320', '21', '2916', '1187', '1482', '2515', '558', '709', '2523', '3438', '2353', '2007', '584', '2814', '2973', '345', '1765', '2040', '837', '1324', '2388', '3416', '3191', '2630', '1870', '788', '1633', '733', '2620', '3623', '2932', '1304', '3194', '318', '740', '1342', '2417', '2350', '598', '654', '3354', '3232', '3755', '2513', '195', '3657', '3341', '2857', '2917', '583', '2640', '3799', '2011', '1735', '3262', '3391', '3056', '3650', '587', '1326', '3196', '175', '2911', '4', '2', '1833', '2491', '1454', '1088', '2647', '1078', '2785', '387', '870', '7', '1816', '3506', '3671', '1923', '2687', '1494', '1248', '3624', '3798', '2686', '3034', '3576', '1565', '1726', '2853', '1155', '1022', '1093', '3340', '1663', '644', '1069', '3033', '1377', '995', '3421', '1284', '1079', '2714', '1706', '1016', '1151', '3030', '695', '1380', '1056', '1826', '739', '2255', '1957', '1430', '1258', '2277', '2147', '1310', '1196', '1658', '1499', '1101', '2724', '2610', '3658', '1013', '622', '257', '1439', '436', '2991', '663', '3692', '1472', '3384', '1405', '658', '2920', '3578', '2615', '2258', '1198', '1808', '1904', '1065', '1285', '2386', '3468', '668', '910', '2493', '1539', '80', '3311', '3660', '3793', '3786', '274', '3118', '3179', '2984', '2397', '2532', '3507', '2614', '2348', '3693', '2546', '2028', '3414', '3165', '738', '277', '2236', '2324', '3697', '763', '720', '57', '2165', '2247', '830', '965', '571', '2167', '2223', '8', '2404', '878', '2235', '1249', '2835', '3383', '749', '3437', '877', '3498', '401', '3412', '2370', '1066', '887', '1037', '1055', '3328', '3379', '265', '409', '3142', '2419', '2918', '2938', '2863', '1156', '1767', '730', '2979', '220', '1974', '2271', '3393', '70', '1961', '578', '1790', '2447', '187', '604', '888', '1766', '3242', '3599', '814', '2650', '1585', '1533', '2125', '3695', '2455', '2466', '3656', '2129', '873', '685', '1321', '2872', '979', '3388', '944', '2854', '2413', '392', '724', '3767', '1921', '1685', '1869', '2234', '388', '2422', '3724', '2293', '1280', '2144', '340', '3590', '3593', '2937', '2598', '298', '2373', '1293', '2203', '3324', '1408', '1948', '2420', '319', '771', '1387', '1091', '3783', '765', '3268', '2140', '1252', '3332', '313', '3069', '1797', '651', '193', '3367', '1087', '1576', '3441', '85', '1357', '2894', '563', '2995', '2209', '3652', '3730', '3404', '858', '2377', '1733', '381', '1235', '2862', '405', '3119', '606', '3329', '2068', '3190', '898', '3433', '1265', '1809', '2059', '2482', '2182', '2711', '2506', '3713', '1389', '2319', '1691', '3223', '3679', '3615', '207', '2763', '2492', '1104', '1532', '824', '412', '1207', '3181', '430', '307', '827', '1577', '2232', '2214', '2196', '3748', '389', '3222', '3235', '2021', '1152', '732', '3154', '1599', '600', '3250', '2569', '967', '2815', '727', '954', '18', '2864', '1323', '2246', '869', '1075', '1234', '2975', '2972', '3435', '3323', '3398', '3325', '1742', '2013', '87', '3319', '2345', '2362', '2505', '3420', '1135', '2721', '2861', '2398', '3675', '2522', '1202', '1144', '3440', '3581', '1984', '2188', '1360', '3318', '3411', '812', '2744', '700', '1891', '1154', '1708', '689', '192', '2208', '3003', '1500', '3770', '2410', '1001', '367', '538', '2900', '1057', '871', '2621', '2475', '1629', '1857', '980', '589', '554', '1206', '2141', '696', '324', '2131', '2990', '3490', '1527', '82', '1329', '364', '1827', '3757', '3442', '2960', '2921', '374', '1774', '1420', '1769', '637', '2789', '3654', '1953', '3145', '2233', '3640', '773', '2923', '3155', '2574', '2070', '3733', '1932', '1024', '840', '796', '1281', '2400', '1391', '3772', '716', '25', '3622', '2860', '189', '1021', '1344', '1388', '2502', '1866', '2446', '3430', '2986', '3041', '2222', '16', '613', '609', '2838', '1873', '2257', '1190', '2197', '73', '413', '835', '3252', '3768', '3719', '1924', '1466', '2531', '1962', '2524', '953', '1806', '876', '3659', '1291', '1730', '581', '1270', '391', '71', '623', '2178', '9', '3394', '424', '3484', '3432', '3765', '1398', '2212', '925', '3510', '1226', '754', '1707', '1739', '1325', '2985', '1764', '842', '2067', '839', '3025', '886', '912', '2479', '406', '2193', '170', '3672', '872', '336', '2406', '1002', '2836', '1197', '524', '2363', '420', '2347', '2393', '1092', '3471', '462', '1945', '2632', '356', '2391', '1910', '2600', '1794', '2858', '2395', '2939', '3415', '1259', '390', '2927', '3620', '1073', '3320', '3618', '3422', '1419', '2148', '3192', '212', '1985', '747', '1530', '3673', '683', '676', '1803', '320', '1282', '1947', '1727', '3603', '225', '1798', '3592', '3579', '479', '84', '2870', '62', '384', '2411', '1095', '1359', '3277', '3233', '792', '3246', '3788', '1114', '657', '3149', '3166', '1203', '232', '1907', '1678', '2468', '2839', '3161', '3010', '2871', '279', '3012', '2581', '2897', '3447', '811', '1307', '2992', '861', '659', '3040', '2305', '2206', '1639', '2084', '1268', '3302', '1390', '2033', '216', '299', '2612', '2250', '115', '349', '1531', '2536', '1729', '2207', '2480', '2308', '2875', '2268', '2401', '2046', '3080', '3632', '3710', '1113', '756', '3306', '2430', '1394', '1385', '3587', '59', '519', '879', '2684', '3193', '1528', '2628', '3313', '2855', '3305', '2253', '2361', '1573', '3202', '3751', '1949', '1058', '3060', '2537', '2622', '2798', '2326', '1470', '3787', '3029', '889', '1784', '5', '625', '567', '1306', '301', '1692', '1563', '2176', '1465', '3472', '3266', '1393', '3651', '1410', '2827', '2251', '354', '1836', '981', '742', '1493', '2168', '969', '2490', '3330', '664', '3399', '723', '3508', '1686', '947', '1100', '977', '3156', '221', '3171', '363', '243', '1813', '2229', '565', '1251', '3244', '866', '542', '1020', '3738', '447', '355', '790', '2922', '3200', '3745', '3214', '1594', '3494', '527', '1222', '3722', '2127', '166', '1521', '255', '1908', '2074', '23', '2914', '2069', '400', '846', '3273', '3706', '1552', '1675', '3641', '2291', '60', '845', '2043', '2462', '632', '968', '172', '1090', '22', '1695', '3322', '1788', '2558', '3631', '310', '2132', '2977', '1649', '841', '1455', '1772', '706', '3013', '1821', '867', '1302', '3159', '2467', '675', '1881', '1299', '1651', '3756', '3063', '1443']

	primers = pd.read_csv('primers.csv', sep=';', dtype={'id': str})

	print('Lendo arquivo BLAST')
	df = pd.read_csv('/media/bioinfo/6tb_hdd/03_ELLEN/02_data/primers_references_blast_alignment/blast_output.tsv', sep='\t', names=['qseqid' ,'pident' ,'gaps' ,'qcovs' ,'mismatch' ,'evalue' ,'score' ,'bitscore' ,'sgi' ,'sacc' ,'slen' ,'qstart' ,'qend' ,'sstart' ,'send' ,'qseq' ,'sseq' ,'length' ,'sstrand' ,'frames' ,'stitle' ,'sscinames'], dtype={'qseqid': str})
	# df_test = df.head(200000)
	# df_test.to_csv('blast_test.tsv', sep='\t', index=False)


	# df = pd.read_csv('blast_test.tsv', sep='\t', dtype={'qseqid': str})
	# df[['Primer', 'SeqID']] = df['qseqid'].apply(split_and_duplicate_primer, args='_')

	print('Separando primers e posições')
	df[['Primer', 'SeqID']] = df['qseqid'].str.split('_', n=1, expand=True)
	print('Separando orientações')
	df[['tmp', 'orientation']] = df['frames'].str.split('/', n=1, expand=True)

	print('Obtendo dataframes dos primers')
	primers_alignment_data = get_primers_alignment_data(df['Primer'].unique(), df)
	print('Limpando dataframes dos primers')
	primers_alignment_data = clean_primers_alignment_data(primers_alignment_data, primers)
	
	print('Obtendo pares de primers')
	accepted_pairs = get_primers_pairs(primers_alignment_data)
	print('Salvando pares de primers')
	write_pairs_file(accepted_pairs, primers)

# time blastn -db /media/bioinfo/6tb_hdd/04_Blast_Databases/PrimersReferences/primers_references -query primers.fasta -outfmt "6 qseqid pident qcovs mismatch sgi sacc slen qstart qend sstart send qseq sseq length sstrand frames stitle sscinames" -out /media/bioinfo/6tb_hdd/03_ELLEN/02_data/primers_references_blast_alignment/blast_output2_less_cols.tsv -num_threads 10 -max_target_seqs 12000 -task blastn-short