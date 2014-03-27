
require 'bio-commandeer'
TEST_DATA_DIR = File.join(File.dirname(__FILE__),'data')

describe 'finishm visualise' do
  path_to_script = File.join(File.dirname(__FILE__),'..','bin','finishm visualise')

  it 'should visualise an entire graph' do
    command = "#{path_to_script} --quiet --fasta #{TEST_DATA_DIR}/explore/1/2seqs.sammy.fa --assembly-dot /dev/stdout"
    stdout = Bio::Commandeer.run(command)
    stdout.should == <<END_OF_DOT
digraph G {
	graph [overlap=scale];
	node [label="\\N"];
	graph [bb="0,0,1533,556"];
	1 [label=n1_length527_coverage21, pos="142,389", width="3.9444", height="0.51389"];
	2 [label=n2_length264_coverage33, pos="142,463", width="3.9444", height="0.51389"];
	3 [label=n3_length227_coverage62, pos="397,537", width="3.9444", height="0.51389"];
	4 [label=n4_length169_coverage16, pos="1077,167", width="3.9444", height="0.51389"];
	5 [label=n5_length155_coverage32, pos="482,463", width="3.9444", height="0.51389"];
	6 [label=n6_length11_coverage71, pos="549,389", width="3.7778", height="0.51389"];
	7 [label=n7_length87_coverage29, pos="781,167", width="3.7778", height="0.51389"];
	8 [label=n8_length102_coverage16, pos="1382,315", width="3.9444", height="0.51389"];
	9 [label=n9_length117_coverage36, pos="1231,241", width="3.9444", height="0.51389"];
	10 [label=n10_length89_coverage71, pos="617,241", width="3.9444", height="0.51389"];
	11 [label=n11_length103_coverage34, pos="441,167", width="4.1111", height="0.51389"];
	12 [label=n12_length11_coverage41, pos="932,93", width="3.9444", height="0.51389"];
	13 [label=n13_length155_coverage38, pos="790,463", width="4.1111", height="0.51389"];
	14 [label=n14_length169_coverage20, pos="1385,167", width="4.1111", height="0.51389"];
	15 [label=n15_length82_coverage21, pos="1080,315", width="3.9444", height="0.51389"];
	16 [label=n16_length103_coverage35, pos="617,315", width="4.1111", height="0.51389"];
	17 [label=n17_length264_coverage19, pos="224,315", width="4.1111", height="0.51389"];
	18 [label=n18_length87_coverage36, pos="774,19", width="3.9444", height="0.51389"];
	2 -> 1 [pos="e,142,407.67 142,444.33 142,436.26 142,426.65 142,417.71"];
	1 -> 17 [color=blue, pos="e,203.6,333.41 162.69,370.33 172.77,361.23 185.04,350.16 195.95,340.32"];
	3 -> 2 [pos="e,200.37,479.94 338.82,520.12 300.52,509 250.24,494.41 210.09,482.76"];
	3 -> 17 [arrowhead=none, pos="368.49,518.69 355.51,509.11 340.88,496.39 331,482 301.25,438.66 323.22,413.01 293,370 282.78,355.46 267.66,342.8 254.17,333.31"];
	3 -> 5 [color=blue, pos="e,460.85,481.41 418.45,518.33 429,509.14 441.86,497.94 453.25,488.03"];
	3 -> 13 [color=blue, pos="e,708.14,478.41 477.95,521.76 541.96,509.7 631.35,492.87 698.05,480.31"];
	4 -> 12 [pos="e,967.36,111.04 1041.5,148.9 1021.9,138.87 997.28,126.31 976.37,115.64"];
	9 -> 4 [pos="e,1114.2,184.88 1193.7,223.09 1172.5,212.9 1145.8,200.05 1123.3,189.22"];
	5 -> 6 [pos="e,532.33,407.41 498.91,444.33 506.98,435.41 516.78,424.59 525.56,414.89"];
	6 -> 11 [arrowhead=none, pos="504.23,371.47 487.67,362.71 470.41,350.41 460,334 430.38,287.31 434.2,218.23 438.18,185.5"];
	6 -> 16 [pos="e,600.08,333.41 566.16,370.33 574.36,361.41 584.3,350.59 593.21,340.89"];
	13 -> 6 [pos="e,604.18,405.94 734.1,445.84 698.24,434.82 651.53,420.48 614,408.96"];
	7 -> 12 [arrowhead=none, pos="817.55,149.09 841.12,137.54 871.66,122.57 895.26,111"];
	10 -> 7 [pos="e,741.67,184.75 656.7,223.09 679.46,212.82 708.22,199.84 732.33,188.96"];
	8 -> 9 [pos="e,1267.7,259 1345.4,297.09 1324.8,286.98 1298.9,274.26 1276.9,263.49"];
	9 -> 14 [pos="e,1347.5,185 1268.3,223.09 1289.4,212.94 1316,200.15 1338.5,189.34"];
	15 -> 9 [pos="e,1194.3,259 1116.6,297.09 1137.2,286.98 1163.1,274.26 1185.1,263.49"];
	10 -> 18 [arrowhead=none, pos="617.51,222.38 618.76,202.73 622.88,171.44 636,148 662.91,99.923 713.93,59.391 745.83,37.207"];
	16 -> 10 [pos="e,617,259.67 617,296.33 617,288.26 617,278.65 617,269.71"];
	10 -> 11 [color=blue, pos="e,483.42,184.84 574.84,223.28 550.21,212.92 518.92,199.76 492.79,188.78"];
	12 -> 18 [pos="e,812.17,36.876 893.75,75.087 872,64.898 844.57,52.052 821.45,41.224"];
	14 -> 12 [pos="e,1020.7,107.49 1295,152.29 1218.6,139.82 1109.3,121.97 1030.6,109.11"];
}
END_OF_DOT
  end

  it 'should visualise particular nodes of interest' do
    command = "#{path_to_script} --quiet --fasta #{TEST_DATA_DIR}/explore/1/2seqs.sammy.fa --assembly-dot /dev/stdout --node-ids 1,17 --leash-length 300"
    stdout = Bio::Commandeer.run(command)
    stdout.should == <<END_OF_DOT
digraph G {
	graph [overlap=scale];
	node [label="\\N"];
	graph [bb="0,0,375,260"];
	1 [label=n1_length527_coverage21, color=red, pos="142,93", width="3.9444", height="0.51389"];
	2 [label=n2_length264_coverage33, pos="142,167", width="3.9444", height="0.51389"];
	3 [label=n3_length227_coverage62, pos="227,241", width="3.9444", height="0.51389"];
	17 [label=n17_length264_coverage19, color=red, pos="227,19", width="4.1111", height="0.51389"];
	2 -> 1 [pos="e,142,111.67 142,148.33 142,140.26 142,130.65 142,121.71"];
	1 -> 17 [color=blue, pos="e,205.85,37.411 163.45,74.327 174,65.143 186.86,53.943 198.25,44.029"];
	3 -> 2 [pos="e,163.15,185.41 205.55,222.33 195,213.14 182.14,201.94 170.75,192.03"];
	3 -> 17 [arrowhead=none, pos="259.08,222.98 272.15,213.77 285.85,201.26 293,186 314.11,140.92 314.11,119.08 293,74 285.91,58.872 272.4,46.454 259.43,37.27"];
}
END_OF_DOT
  end
end
