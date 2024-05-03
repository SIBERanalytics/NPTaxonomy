# Import necessary packages
import pandas as pd
import json
import matplotlib.pyplot as plt 
plt.rcParams['font.family'] = ['Arial']
plt.rcParams.update({'font.size':16})

# Organize t-SNE results
results_map4 = []
results_mpn = []
results_lastFFN = []
for i in range(4):
    # Fingerprint MAP4
    with open(f'MAP4_1E{i+2}.json', 'r') as f:
        results_map4_i = json.load(f)
    results_map4.append(results_map4_i)
    # Fingerprint MPN   
    with open(f'MPN_1E{i+2}.json', 'r') as f:
        results_mpn_i = json.load(f)       
    results_mpn.append(results_mpn_i)   
    # Fingerprint last_FFN
    with open(f'last_FFN_1E{i+2}.json', 'r') as f:
        results_lastFFN_i = json.load(f)       
    results_lastFFN.append(results_lastFFN_i) 
# Combine all dictionaries to dataframes
results_map4 = pd.DataFrame(results_map4)
results_mpn = pd.DataFrame(results_mpn)
results_lastFFN = pd.DataFrame(results_lastFFN)
                  
# Figure of divergence against perplexity
plt.figure()
plt.semilogx(results_map4['tsne_perplexity'],results_map4['tsne_kl_divergence'], 'g-*',label='MAP4')
plt.semilogx(results_mpn['tsne_perplexity'],results_mpn['tsne_kl_divergence'], 'r-*',label='MPN')
plt.semilogx(results_lastFFN['tsne_perplexity'],results_lastFFN['tsne_kl_divergence'], 'b-*',label='last_FFN')
plt.ylabel('KL divergence'),plt.xlabel('Perplexity')
plt.ylim([0,3])
plt.legend()
plt.tight_layout()
plt.savefig('kl_divergence.jpg', dpi=300)

# Figure of db_score against perplexity
plt.figure()
plt.semilogx(results_map4['tsne_perplexity'],results_map4['db_score'], 'g-*',label='MAP4')
plt.semilogx(results_mpn['tsne_perplexity'],results_mpn['db_score'], 'r-*',label='MPN')
plt.semilogx(results_lastFFN['tsne_perplexity'],results_lastFFN['db_score'], 'b-*',label='last_FFN')
plt.ylabel('DB score'),plt.xlabel('Perplexity')
plt.legend()
plt.tight_layout()
plt.savefig('db_score.jpg', dpi=300)



