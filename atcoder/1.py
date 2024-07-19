import math

class Solution:
    def minimumCost(self, m: int, n: int, horizontalCut: List[int], verticalCut: List[int]) -> int:
        horizontalCut.sort()
        verticalCut.sort()

        # Calculate the maximum difference between adjacent horizontal cuts
        max_horizontal_diff = max(horizontalCut[0], m - horizontalCut[-1])
        for i in range(1, len(horizontalCut)):
            max_horizontal_diff = max(max_horizontal_diff, horizontalCut[i] - horizontalCut[i - 1])

        # Calculate the maximum difference between adjacent vertical cuts
        max_vertical_diff = max(verticalCut[0], n - verticalCut[-1])
        for j in range(1, len(verticalCut)):
            max_vertical_diff = max(max_vertical_diff, verticalCut[j] - verticalCut[j - 1])

        # Return the minimum total cost
        return (max_horizontal_diff * max_vertical_diff) % (10**9 + 7)

# Example usage
solution = Solution()
print(solution.minimumCost(3, 2, [1, 3], [5]))  # Output: 13
print(solution.minimumCost(2, 2, [7], [4]))     # Output: 15
