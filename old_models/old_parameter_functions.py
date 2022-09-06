'''These functions have been replaced with newer version in parameters.py or have become obsolete'''

def generate_loo_interaction_matrices(interactions):
    n = len(interactions[0])
    loo_interactions = []
    for i in range(0, n):
        temp_interactions = np.delete(interactions, i, 0)
        temp_interactions = np.delete(temp_interactions, i, 1)
        loo_interactions.append(temp_interactions)
    return np.array(loo_interactions)