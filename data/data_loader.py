from ucimlrepo import fetch_ucirepo
import pandas as pd
import random

def load_data_adults():
    """Lengths and preprocesses the dataset."""
    print("Loading dataset...")
    data = fetch_ucirepo(id=2)
    data = data["data"]["features"][
        [
            "age",
            "education-num",
            "sex",
            "capital-gain",
            "capital-loss",
            "hours-per-week",
        ]
    ]
    map_sex = {"Male": 0, "Female": 1}
    data["sex"] = data["sex"].map(map_sex)

    # Drop duplicate rows

    normalized_data = data.apply(lambda x: (x - x.min()) / (x.max() - x.min()))
    normalized_data = normalized_data.drop_duplicates()
    print(f"Data loaded. Total records: {normalized_data.shape[0]}")
    
    # split into v and u
    female = normalized_data[normalized_data["sex"] == 1]
    male = normalized_data[normalized_data["sex"] == 0]
    
    female = female.drop("sex", axis=1).drop_duplicates()
    male = male.drop("sex", axis=1).drop_duplicates()
    
    merged = male.merge(female, how='left', indicator=True)

    male = merged[merged['_merge'] == 'left_only'].drop(columns='_merge')
    
    v = [tuple(i) for i in female.to_numpy().tolist()]
    u = [tuple(i) for i in male.to_numpy().tolist()]
    
    return normalized_data, v, u

def load_data_credits():
    """Lengths and preprocesses the dataset."""
    print("Loading dataset...")
    data = fetch_ucirepo(id=350)
    data = data["data"]["features"][
        [
            "X2",
            "X12",
            "X13",
            "X14",
            "X15",
            "X16",
            "X17",
        ]
    ]
    cols_normalize = [ "X12", "X13", "X14", "X15", "X16", "X17"]

    for col in cols_normalize:
        data[col] = (data[col] - data[col].min()) / (data[col].max() - data[col].min())

    normalized_data = data.drop_duplicates()
    print(f"Data loaded. Total records: {normalized_data.shape[0]}")
    
    # split into v and u
    female = normalized_data[normalized_data["X2"] == 2]
    male = normalized_data[normalized_data["X2"] == 1]

    female = female.drop("X2", axis=1).drop_duplicates()
    male = male.drop("X2", axis=1).drop_duplicates()
    merged = male.merge(female, how='left', indicator=True)

    male = merged[merged['_merge'] == 'left_only'].drop(columns='_merge')

    v = [tuple(i) for i in male.to_numpy().tolist()]
    u = [tuple(i) for i in female.to_numpy().tolist()]
    
    return normalized_data, v, u

def load_data_gamma():
    """Lengths and preprocesses the dataset."""
    print("Loading dataset...")
    data = fetch_ucirepo(id=159)
    data = data["data"]["original"]
    map_class = {"g": 0, "h": 1}
    data["class"] = data["class"].map(map_class)

    normalized_data = data.apply(lambda x: (x - x.min()) / (x.max() - x.min()))
    normalized_data = normalized_data.drop_duplicates()
    print(f"Data loaded. Total records: {normalized_data.shape[0]}")
    
    # split into v and u
    female = normalized_data[normalized_data["class"] == 1]
    male = normalized_data[normalized_data["class"] == 0]

    female = female.drop("class", axis=1).drop_duplicates()
    male = male.drop("class", axis=1).drop_duplicates()
    merged = male.merge(female, how="left", indicator=True)

    male = merged[merged["_merge"] == "left_only"].drop(columns="_merge")

    u = [tuple(i) for i in male.to_numpy().tolist()]
    v = [tuple(i) for i in female.to_numpy().tolist()]
    
    return normalized_data, v, u

def load_data_popsim(seed, v_min, v_max, u_min, u_max):
    """Lengths and preprocesses the dataset."""
    print("Loading dataset...")
    data = pd.read_csv("data/popsim_1M.csv", usecols=["race", "lon", "lat"])

    most_popular_races = data['race'].value_counts().nlargest(2).index
    data = data[data['race'].isin(most_popular_races)]

    data['race'] = data['race'].astype('category').cat.codes

    data = data.astype(float)

    normalized_data = data.apply(lambda x: (x - x.min()) / (x.max() - x.min()))
    normalized_data = normalized_data.drop_duplicates()
    
    # split into v and u
    print(f"Data loaded. Total records: {normalized_data.shape[0]}")
    v_all = normalized_data[normalized_data["race"] == 1]
    u_all = normalized_data[normalized_data["race"] == 0]
    
    v_all = v_all.drop("race", axis=1).drop_duplicates()
    u_all = u_all.drop("race", axis=1).drop_duplicates()
    
    merged = v_all.merge(u_all, how='left', indicator=True)
    v_all = merged[merged['_merge'] == 'left_only'].drop(columns='_merge')
    merged = u_all.merge(v_all, how='left', indicator=True)
    u_all = merged[merged['_merge'] == 'left_only'].drop(columns='_merge')
    
    v_all = [tuple(i) for i in v_all.to_numpy().tolist()]
    u_all = [tuple(i) for i in u_all.to_numpy().tolist()]
    
    random.seed(seed)
    random.shuffle(v_all)
    random.shuffle(u_all)
    
    v_size = random.randint(v_min, v_max)
    u_size = random.randint(u_min, u_max)
    v = v_all[:v_size]
    u = u_all[:u_size]
    
    return normalized_data, v, u