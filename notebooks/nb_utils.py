

def remove_rogue_tqdm():
    import tqdm

    try:
        tqdm._instances.clear()
    except AttributeError:
        pass


def describe(df, random=False):
    print(df.shape)
    print("--- First 5 entries ---")
    display(df.head())
    if random:
        print("--- Random subset ---")
        display(df.sample(5))
