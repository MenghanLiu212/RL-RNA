import argparse
import util.rna_toolkit
import util.pseudoknotted
import util.pseudoknot_free

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--real_rna_path',
        type = str,
        required = True,
        help = "Real RNA path"
    )

    parser.add_argument(
        '--simulation_rna_path',
        type = str,
        required = True,
        help = "Simulation RNA path"
    )

    parser.add_argument(
        '--is_pseudoknot_free',
        type = str,
        required = True,
        help = "is it training for pseduoknot_free RNAs(Y/n)?"
    )
    
    args = parser.parse_args()
    real_rna_path = args.real_rna_path
    simulation_rna_path = args.simulation_rna_path
    pseudoknot_free_check = args.is_pseudoknot_free.lower()

    if pseudoknot_free_check != 'y':
        pseudonkot_message = "Training ENTRNA for pseudoknotted RNA"
        accuracy = util.pseudoknotted.entrna_train(real_rna_path, simulation_rna_path)
    else:
        pseudonkot_message = "Training ENTRNA pseudoknot-free RNA"
        accuracy = util.pseudoknot_free.entrna_train(real_rna_path, simulation_rna_path)
    
    print '\n\n\n\n===============================================================\n\n'
    print pseudonkot_message
    print 'ENTRNA training accuracy:',accuracy
    print '\n\n===============================================================\n\n'
