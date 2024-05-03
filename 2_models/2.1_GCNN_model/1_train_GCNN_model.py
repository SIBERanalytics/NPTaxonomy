# Import necessary package
import chemprop

# Define arguments
arguments = [
    '--data_path', 'data/np_5classes.csv',
    '--separate_val_path','data/np_5classes.csv',   
    
    '--dataset_type', 'multiclass',
    '--multiclass_num_classes', '5',
    '--extra_metrics' , 'accuracy','f1','mcc', 
    '--split_sizes', '1', '0', '0',
    '--epochs', '200',
    
    '--save_dir', 'results',   
    '--config_path', 'results.json',
    '--num_workers', '16'    
]
args = chemprop.args.TrainArgs().parse_args(arguments)

# Train GCNN model
mean_score, std_score  = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)
