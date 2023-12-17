#Eden Ahady 318948106
import pandas as pd


class myData:
    #constarctor
    def __init__(self, path) -> None:
        self.data = pd.read_csv(path)

    def num_players_height(self, x, y):
        """

        :param x: is one height
        :param y: i second height
        :return: the numbers of player their height is between x and y
        """
        #list of heights
        heights = self.data.loc[:, "height_cm"]
        count = 0;
        for i in range(len(heights)):
            if x <= heights[i] <= y:
                count = count + 1
        return(count)

    def df_birthyear(self, year):
        """

        :param year: is the year given
        :return: data frame of the players who were born in this year
        """
        # convert the "date_of_birth" column to a datetime data type
        self.data["dob"] = pd.to_datetime(self.data["dob"])

        # boolean mask that selects where the year of the "date_of_birth" column is 1980
        mask = self.data["dob"].dt.year == year

        # rows where the year of the "date_of_birth" column is 1980,"short_name" and "club" columns
        df_1980 = self.data.loc[mask, ["short_name", "club_name"]]

        return(df_1980)

    def list_sorted(self, name1, name2, k):
        """

        :param name1: column mane
        :param name2: column mane
        :param k: number of players we want to return
        :return: list of short name of the k highest players, sorted by 'name2'
        """
        #first highest k values
        filter1 = self.data.sort_values(name1, ascending=False).head(k)
        #sort by the second column given
        filter2 = filter1.sort_values(name2, ascending=True)
        return(filter2["short_name"].tolist())

    def tuples_players_by_year(self, year1, year2):
        """

        :param year1: is year given
        :param year2: is year given
        :return: list of tuples which contains the year and the number of players who born this year
        """
        self.data["dob"] = pd.to_datetime(self.data["dob"])

        # group the DataFrame by "date_of_birth" column
        year_group = self.data.groupby(self.data["dob"].dt.year)

        # list of tuples that holds the year and the number of players born in that year
        players_and_year = [(year, len(group)) for year, group in year_group]

        # only years between year_start and year_end
        filter_list = [x for x in players_and_year if x[0] <= year2 and x[0] >= year1]
        return(filter_list)

    def mean_std(self, column, name):
        """

        :param column: name of column
        :param name: first name
        :return:mean and std of the values in the column for those players
        """
        mean_name = name + " "
        #getting the column that starts with the name given
        mean_sort = self.data["long_name"].str.match(mean_name)
        our_data = self.data[mean_sort]
        # calculate the mean and standard deviation of the values in the given column
        mean = our_data[column].mean()
        std = our_data[column].std()
        return(mean, std)

    def max_players(self, column):
        """

        :param column: name of column
        :return:the value which contains the maximum players in this column
        """
        #getting the value from the column with the maximum players
        maximum = self.data[column].value_counts().idxmax()
        return(maximum)




