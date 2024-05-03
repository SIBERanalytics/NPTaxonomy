# Import necessary package
import chemprop

# MPN fingerprints of 133092 NPs from LOTUS database
arguments1 = [
    '--test_path', 'data/np_5classes.csv',
    '--preds_path', 'lotus_MPN.csv',
    '--checkpoint_dir', 'results',
    '--num_workers', '16',
    '--fingerprint_type', 'MPN']
    
args1 = chemprop.args.FingerprintArgs().parse_args(arguments1)
lotus_MPN  = chemprop.train.molecule_fingerprint.molecule_fingerprint(args=args1)

# last_FFN fingerprints of 133092 NPs from LOTUS database
arguments2 = [
    '--test_path', 'data/np_5classes.csv',
    '--preds_path', 'lotus_last_FFN.csv',
    '--checkpoint_dir', 'results',
    '--num_workers', '16',
    '--fingerprint_type', 'last_FFN']
    
args2 = chemprop.args.FingerprintArgs().parse_args(arguments2)
lotus_last_FFN = chemprop.train.molecule_fingerprint.molecule_fingerprint(args=args2)
