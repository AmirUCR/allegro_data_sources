import os
import sys
from downloaders.NCBI.ncbi_downloader import NCBI_Downloader
from downloaders.FungiDB.fungidb_downloader import FungiDB_Downloader
from downloaders.EnsemblFungi.ensembl_download import EnsemblFungi_Downloader
from downloaders.MycoCosm.mycocosm_download import MycoCosm_Downloader
from utils.merger import merge_dbs


def main(choice_arg=''):
    if not os.path.exists('data'):
        os.mkdir('data')

    choice = '0'

    downloaders = {
        '1': NCBI_Downloader,
        '2': FungiDB_Downloader,
        '3': EnsemblFungi_Downloader,
        '4': MycoCosm_Downloader,
    }

    if choice_arg:
        print(f'Received choice {choice_arg} via argv.')
        choice = choice_arg
    else:
        print('Welcome to the ALLEGRO dataset downloader. Select a download option:\n')
        choice = input('1. NCBI Datasets\n2. FungiDB\n3. EnsemblFungi\n4. MycoCosm\n5. Merge Downloads\n6. Quit\nYour choice: ')

    if choice in downloaders:
        downloader = downloaders[choice]()
        downloader.download()
    elif choice == '6':
        merge_dbs()

    print('Goodbye.')
    return 0
    

if __name__ == '__main__':
    if len(sys.argv) > 1:
        sys.exit(main(sys.argv[1]))
    sys.exit(main())