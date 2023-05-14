# Hold results of calculations. (Used for plotting results)
class Holder:
    def __init__(self, time_array, response_array, response_edited, response_edited_with_k, af_array):
        self.time_array = time_array
        self.response_array = response_array
        self.af_array = af_array
        self.response_edited = response_edited
        self.response_edited_with_k = response_edited_with_k