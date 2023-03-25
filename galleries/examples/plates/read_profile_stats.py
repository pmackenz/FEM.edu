import pstats
from pstats import SortKey

p = pstats.Stats('profile_data_full.txt')
p.strip_dirs() #.sort_stats(-1).print_stats()
p.sort_stats(SortKey.NAME)
#p.print_stats()

p.sort_stats(SortKey.CUMULATIVE).print_stats(20)

p.sort_stats(SortKey.TIME).print_stats(20)


p = pstats.Stats('profile_data_sparse.txt')
p.strip_dirs() #.sort_stats(-1).print_stats()
p.sort_stats(SortKey.NAME)
#p.print_stats()

p.sort_stats(SortKey.CUMULATIVE).print_stats(20)

p.sort_stats(SortKey.TIME).print_stats(20)
