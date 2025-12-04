import sys
import os
sys.path.append(os.getcwd())
from tf_explorer import encode_client

print("Searching for CTCF...")
results = encode_client.search_encode_tf_chipseq('CTCF', 'Homo sapiens')
print(f"Found {len(results)} results.")
if results:
    print("First result:", results[0])
else:
    print("No results found.")
