import argparse
import util.rna_toolkit
import util.pseudoknotted
import util.pseudoknot_free

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--seq_file',
        type = str,
        required = True,
        help = "RNA sequence file path"
    )

    parser.add_argument(
        '--str_file',
        type = str,
        required = True,
        help = "RNA secondary structure file path"
    )
    
    args = parser.parse_args()
    rna_seq = open(args.seq_file).read()
    rna_str = open(args.str_file).read()[1:-1].split(",")
    rna_str = [int(i) for i in rna_str]

    if util.rna_toolkit.is_pseudoknotted(rna_str) > 0:
        pseudonkot_message = "This is pseudoknotted RNA"
        foldability = util.pseudoknotted.entrna_main(rna_seq, rna_str)
    else:
        pseudonkot_message = "This is pseudoknot-free RNA"
        foldability = util.pseudoknot_free.entrna_main(rna_seq, rna_str)
    
    print '\n\n\n\n===============================================================\n\n'
    print 'RNA sequence:\n',rna_seq
    print 'RNA secondary structure:\n', rna_str
    print pseudonkot_message
    print 'Foldability:',foldability
    print '\n\n===============================================================\n\n'
