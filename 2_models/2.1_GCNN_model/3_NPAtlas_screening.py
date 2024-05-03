# Import necessary package
import chemprop

# GCNN prediction for NPAtlas dataset (13136 NPs)
arguments = [
    '--test_path', 'npatlas_dataset.csv',
    '--preds_path', 'npatlas_preds.csv',
    '--checkpoint_dir', 'results'
]
args = chemprop.args.PredictArgs().parse_args(arguments)
npatlas_preds = chemprop.train.make_predictions(args=args)

# MPN fingerprints of 13136 NPs from NPAtlas dataset
arguments1 = [
    '--test_path', 'npatlas_dataset.csv',
    '--preds_path', 'npatlas_MPN.csv',
    '--checkpoint_dir', 'results',
    '--fingerprint_type', 'MPN']
args1 = chemprop.args.FingerprintArgs().parse_args(arguments1)
npatlas_MPN  = chemprop.train.molecule_fingerprint.molecule_fingerprint(args=args1)

# last_FFN fingerprints of 13136 NPs from NPAtlas dataset
arguments2 = [
    '--test_path', 'npatlas_dataset.csv',
    '--preds_path', 'npatlas_last_FFN.csv',
    '--checkpoint_dir', 'results',
    '--fingerprint_type', 'last_FFN']
args2 = chemprop.args.FingerprintArgs().parse_args(arguments2)
npatlas_last_FFN  = chemprop.train.molecule_fingerprint.molecule_fingerprint(args=args2)
