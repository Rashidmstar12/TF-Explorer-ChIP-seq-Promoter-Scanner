import sys
print("Python version:", sys.version)
try:
    import requests
    print("Requests version:", requests.__version__)
except ImportError:
    print("Requests not installed")

import os
sys.path.append(os.getcwd())

try:
    from tf_explorer import encode_client
    print("Imported encode_client successfully")
    results = encode_client.search_encode_tf_chipseq('CTCF', 'Homo sapiens')
    print("Results count:", len(results))
except Exception as e:
    print("Error:", e)
