#!/usr/bin/env python

import argparse
import requests
import json
import re

def fetch_pfam_description(pfam_id):
    """Fetch the description for a given Pfam ID from the InterPro API."""
    url = f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}'
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.json()
    else:
        print(f'Error retrieving data for {pfam_id}: {response.status_code}')
        return None

def main(pfam_ids):
    """Main function to process Pfam IDs and fetch their descriptions."""
    for pfam_id in pfam_ids:
        # Trim suffix if it exists
        trimmed_id = re.sub(r'\.\d+$', '', pfam_id)
        
        # Fetch description
        data = fetch_pfam_description(trimmed_id)
        
        if data:
            name = data['metadata']['name']
            description = data['metadata'].get('description', 'No description available')
            output = {
                'Pfam ID': trimmed_id,
                'Name': name,
                'Description': description
            }
            print(json.dumps(output, indent=4))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch functional descriptions for Pfam IDs.')
    parser.add_argument('-i', '--ids', nargs='+', required=True, help='List of Pfam IDs (at least one).')
    
    args = parser.parse_args()
    
    main(args.ids)